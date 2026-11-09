// Status TLVs for the time-station tick detector: encode and decode in one place
// Copyright 2026, HamSCI
#include <math.h>
#include <string.h>
#include "status.h"
#include "pilot.h"
#include "tick_status.h"

// scale (float, big endian) then n int16 values, big endian. Scale is the
// largest magnitude, so values use the full int16 range and the quantisation
// error is at most scale/65534. Returns the byte count
static size_t pack_scaled(uint8_t *out,float const *v,int n,float scale){
  uint32_t u;
  memcpy(&u,&scale,sizeof u);
  for(int i=0; i < 4; i++)
    out[i] = (uint8_t)(u >> (24 - 8*i));
  uint8_t *cp = out + 4;
  for(int i=0; i < n; i++){
    long q = scale > 0 ? lrintf(v[i] / scale * 32767.0f) : 0;
    if(q > 32767) q = 32767;
    if(q < -32767) q = -32767;
    *cp++ = (uint8_t)((uint16_t)q >> 8);
    *cp++ = (uint8_t)q;
  }
  return cp - out;
}

static float max_abs(float const *v,int n){
  float m = 0;
  for(int i=0; i < n; i++)
    if(fabsf(v[i]) > m)
      m = fabsf(v[i]);
  return m;
}

void encode_tick_status(uint8_t **bp,struct tick_second const * const s,int64_t const second_utc_ns,
                        bool const used_pilot,unsigned const templates,int64_t const base_input,int const decim){
  encode_int64(bp,TICK_SECOND_UTC,(uint64_t)second_utc_ns);
  encode_int64(bp,TICK_WINDOW_START,(uint64_t)pilot_out_to_input(base_input,s->window_start_index,decim,NULL));
  encode_int(bp,TICK_STEP,s->step * decim);
  encode_int(bp,TICK_FLAGS,(s->complete ? 1 : 0) | (used_pilot ? 2 : 0) | (int)((templates & 7) << 4));
  static enum status_type const ptype[TICK_NTEMPLATES] = {TICK_PROFILE_WWV,TICK_PROFILE_WWVH,TICK_PROFILE_BPM};
  uint8_t buf[4 + 4 * TICK_MAX_PROFILE];
  float v[2 * TICK_MAX_PROFILE];
  float peaks[2 * TICK_NTEMPLATES];
  for(int k=0; k < TICK_NTEMPLATES; k++){
    peaks[k] = NAN;
    peaks[TICK_NTEMPLATES + k] = 0;
    if(!(templates & (1u << k)))
      continue;
    for(int i=0; i < s->nprofile; i++){
      v[2*i] = crealf(s->profile[k][i]);
      v[2*i+1] = cimagf(s->profile[k][i]);
    }
    size_t const n = pack_scaled(buf,v,2 * s->nprofile,max_abs(v,2 * s->nprofile));
    encode_string(bp,ptype[k],buf,n); // opaque bytes; the multi-byte length form covers > 127
    peaks[k] = (float)((s->peak_index[k] - (double)s->window_start_index) * decim); // NaN stays NaN
    peaks[TICK_NTEMPLATES + k] = (float)s->peak_snr_db[k];
  }
  encode_vector(bp,TICK_PEAKS,peaks,2 * TICK_NTEMPLATES);
  float const mk[TICK_MARKER_TONES] = {
    (float)s->marker_energy[0],(float)s->marker_energy[1],(float)s->marker_energy[2] };
  encode_vector(bp,TICK_MARKERS,mk,TICK_MARKER_TONES);
  // The core's int16 envelope goes out as it is, its scale beside it
  uint8_t ebuf[4 + 2 * TICK_ENVELOPE_LEN];
  uint32_t u;
  memcpy(&u,&s->envelope_scale,sizeof u);
  for(int i=0; i < 4; i++)
    ebuf[i] = (uint8_t)(u >> (24 - 8*i));
  for(int i=0; i < TICK_ENVELOPE_LEN; i++){
    ebuf[4 + 2*i] = (uint8_t)((uint16_t)s->envelope[i] >> 8);
    ebuf[5 + 2*i] = (uint8_t)s->envelope[i];
  }
  encode_string(bp,TICK_ENVELOPE,ebuf,sizeof ebuf);
}

static float get_float(uint8_t const *cp){
  uint32_t u = 0;
  for(int i=0; i < 4; i++)
    u = (u << 8) | cp[i];
  float f;
  memcpy(&f,&u,sizeof f);
  return f;
}

static int16_t get_i16(uint8_t const *cp){
  return (int16_t)(uint16_t)((cp[0] << 8) | cp[1]);
}

bool decode_tick_status(struct tick_status * const ts,enum status_type const type,uint8_t const * const cp,int const optlen){
  switch(type){
  case TICK_SECOND_UTC:
    ts->second_utc_ns = (int64_t)decode_int64(cp,optlen);
    ts->valid = true;
    break;
  case TICK_WINDOW_START:
    ts->window_start_input = (int64_t)decode_int64(cp,optlen);
    break;
  case TICK_STEP:
    ts->step_input = decode_int(cp,optlen);
    break;
  case TICK_FLAGS:
    ts->flags = decode_int(cp,optlen);
    ts->complete = ts->flags & 1;
    ts->used_pilot = (ts->flags >> 1) & 1;
    ts->templates = (ts->flags >> 4) & 7;
    break;
  case TICK_PROFILE_WWV:
  case TICK_PROFILE_WWVH:
  case TICK_PROFILE_BPM:
    {
      int const k = type == TICK_PROFILE_WWV ? TICK_WWV : type == TICK_PROFILE_WWVH ? TICK_WWVH : TICK_BPM;
      ts->nprofile[k] = 0;
      if(optlen >= 4 && (optlen - 4) % 4 == 0 && (optlen - 4) / 4 <= TICK_MAX_PROFILE){
        int const n = (optlen - 4) / 4;
        float const scale = get_float(cp);
        ts->profile_scale[k] = scale;
        for(int i=0; i < n; i++)
          ts->profile[k][i] = (get_i16(cp + 4 + 4*i) + I * get_i16(cp + 6 + 4*i)) * (scale / 32767.0f);
        ts->nprofile[k] = n;
      }
    }
    break;
  case TICK_PEAKS:
    if(optlen == 2 * TICK_NTEMPLATES * 4){
      for(int k=0; k < TICK_NTEMPLATES; k++){
        ts->peak_offset[k] = get_float(cp + 4*k);   // may be NaN
        ts->peak_snr_db[k] = get_float(cp + 4*(TICK_NTEMPLATES + k));
      }
      ts->have_peaks = true;
    }
    break;
  case TICK_MARKERS:
    if(optlen == TICK_MARKER_TONES * 4){
      for(int i=0; i < TICK_MARKER_TONES; i++)
        ts->marker_energy[i] = get_float(cp + 4*i);
      ts->have_markers = true;
    }
    break;
  case TICK_ENVELOPE:
    ts->nenvelope = 0;
    if(optlen == 4 + 2 * TICK_ENVELOPE_LEN){
      ts->envelope_scale = get_float(cp);
      for(int i=0; i < TICK_ENVELOPE_LEN; i++)
        ts->envelope[i] = get_i16(cp + 4 + 2*i);
      ts->nenvelope = TICK_ENVELOPE_LEN;
    }
    break;
  default:
    return false;
  }
  return true;
}

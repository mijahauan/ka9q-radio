// Status TLVs for the time-station tick detector: TICK_*
// One encode/decode pair, used by radio_status.c, decode_status.c and the roundtrip test
// Copyright 2026, HamSCI
#ifndef _TICK_STATUS_H
#define _TICK_STATUS_H 1
#include <stdint.h>
#include <stdbool.h>
#include <complex.h>
#include "status.h"
#include "tick.h"

// One report as a packet carries it: positions on the INPUT counter
struct tick_status {
  bool valid;                  // TICK_SECOND_UTC seen
  int64_t second_utc_ns;
  int64_t window_start_input;
  int step_input;              // input samples per profile step
  int flags;
  bool complete, used_pilot;
  unsigned templates;          // mask of 1 << enum tick_template
  int nprofile[TICK_NTEMPLATES];            // 0 = that profile did not arrive
  float profile_scale[TICK_NTEMPLATES];
  float complex profile[TICK_NTEMPLATES][TICK_MAX_PROFILE]; // dequantised
  bool have_peaks;
  float peak_offset[TICK_NTEMPLATES];       // input samples from window start; NaN = none
  float peak_snr_db[TICK_NTEMPLATES];
  bool have_markers;
  float marker_energy[TICK_MARKER_TONES];
  int nenvelope;               // 0 or TICK_ENVELOPE_LEN
  float envelope_scale;
  int16_t envelope[TICK_ENVELOPE_LEN];
};

// Append the TICK_* group for report s. second_utc_ns and used_pilot come from
// the glue's tick_meta, templates is the channel's mask. s's window start is
// on the channel-output counter, (input - base_input) / decim, so the input
// index is pilot_out_to_input(base_input, window_start_index, decim). An
// incomplete report goes out too, with zeros, NaN peaks and flag bit0 clear
void encode_tick_status(uint8_t **bp, struct tick_second const *s, int64_t second_utc_ns,
                        bool used_pilot, unsigned templates, int64_t base_input, int decim);
// Decode one TLV into *ts (zero it before the first). True when type was one of ours
bool decode_tick_status(struct tick_status *ts, enum status_type type, uint8_t const *cp, int optlen);

// Once-only rule: true the first time report number last_seq (counted from 1; 0 = no
// report) is asked about for the destination whose last-sent number is *sent_seq, and
// records it; false for the same report afterwards
static inline bool tick_status_take(uint64_t last_seq,uint64_t *sent_seq){
  if(last_seq == 0 || last_seq == *sent_seq)
    return false;
  *sent_seq = last_seq;
  return true;
}

#endif

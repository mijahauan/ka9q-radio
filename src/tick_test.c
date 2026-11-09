// Unit tests for the time-station tick detector (tick.c)
// Each test synthesizes a few seconds of 24 kHz carrier-centred AM baseband
// and feeds it in 480-sample (20 ms) batches
// make tick_test && ./tick_test ; exits 0 when every test passes
#undef NDEBUG
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tick.h"

#define SR 24000
#define BATCH 480
#define S0 (SR/2)          // first whole second sits half a second into the stream
#define MAXREP 8

// Uniform (0,1) from xorshift64*, then Box-Muller
static uint64_t Rng = 88172645463325252ULL;
static double uniform(void){
  Rng ^= Rng >> 12; Rng ^= Rng << 25; Rng ^= Rng >> 27;
  return ((Rng * 2685821657736338717ULL) >> 11) * 0x1.0p-53 + 0x1.0p-54;
}
static double gauss(void){
  return sqrt(-2*log(uniform())) * cos(2*M_PI*uniform());
}

// One carrier (amplitude 1) AM-modulated by the ticks of up to two stations,
// a marker tone at second 0 of the minute and a steady test tone
struct synth {
  double tick_ms[2], tick_amp[2];   // WWV (1000 Hz, 5 cycles), WWVH (1200 Hz, 6 cycles); amp 0 = absent
  double marker_hz, marker_amp;     // 800 ms at second 0
  double tone_hz, tone_depth;       // steady AM tone
  double noise_rms;                 // complex rms
  int first_som;                    // second-of-minute of the first whole second
  double phase;                     // carrier phase, radians
};

static float complex synth_sample(struct synth const *sy,int64_t n){
  int64_t const rel = n - S0;
  int64_t sec = rel >= 0 ? rel / SR : -((-rel + SR - 1) / SR);
  int64_t const pos = rel - sec * SR;
  int const som = (int)(((sec + sy->first_som) % 60 + 60) % 60);
  double mod = 0;
  static double const freq[2] = {1000,1200};
  static double const cycles[2] = {5,6};
  if(som != 29 && som != 59 && !(som == 0 && sy->marker_amp > 0)){
    for(int k=0; k < 2; k++){
      if(sy->tick_amp[k] == 0)
        continue;
      double const t = (pos - sy->tick_ms[k] * 1e-3 * SR) / SR;
      if(t >= 0 && t < cycles[k] / freq[k])
        mod += sy->tick_amp[k] * cos(2*M_PI*freq[k]*t);
    }
  }
  if(som == 0 && sy->marker_amp > 0 && pos < 0.8 * SR)
    mod += sy->marker_amp * cos(2*M_PI*sy->marker_hz*pos/SR);
  if(sy->tone_depth > 0)
    mod += sy->tone_depth * cos(2*M_PI*sy->tone_hz*n/SR);
  double complex s = cexp(I*sy->phase) * (1 + mod);
  if(sy->noise_rms > 0)
    s += sy->noise_rms * M_SQRT1_2 * (gauss() + I*gauss());
  return (float complex)s;
}

// The nominal second a block belongs to: the one holding its first sample
static int64_t second_of(int64_t index){
  if(index < S0)
    return S0;
  return S0 + ((index - S0) / SR) * SR;
}

static struct tick_second Rep[MAXREP];
static int Ret[MAXREP];

// Feed until `seconds` whole seconds from S0 have reported (plus the half
// second before S0); skip_block >= 0 leaves that block out.
// Returns the number of reports; Ret[] holds the gap returns seen per report
static int run(struct synth const *sy,unsigned mask,int seconds,int skip_block,int *gaps){
  struct tick_config const cfg = { .samprate = SR, .templates = mask };
  struct tick_state *t = tick_create(&cfg);
  assert(t != NULL);
  static float complex buf[BATCH];
  int nrep = 0;
  int64_t const total = S0 + (int64_t)(seconds-1) * SR + tick_span_end(t) + BATCH;
  if(gaps)
    *gaps = 0;
  int block = 0;
  for(int64_t index = 0; index < total; index += BATCH, block++){
    for(int i=0; i < BATCH; i++)
      buf[i] = synth_sample(sy,index + i);
    if(block == skip_block)
      continue;
    struct tick_second r;
    int const ret = tick_process(t,buf,BATCH,index,second_of(index),&r);
    if(ret < 0 && gaps)
      (*gaps)++;
    if(ret == 1 && nrep < MAXREP){
      Rep[nrep] = r;
      Ret[nrep] = ret;
      nrep++;
    }
  }
  tick_destroy(t);
  return nrep;
}

static double maxmag(float complex const *p,int n){
  double m = 0;
  for(int k=0; k < n; k++)
    if(cabsf(p[k]) > m)
      m = cabsf(p[k]);
  return m;
}

static void test_two_stations_separate(void){
  struct synth sy = { .tick_ms = {4.2,29.0}, .tick_amp = {1.0,1.0}, .noise_rms = 0.02,
                      .first_som = 10, .phase = 0.7 };
  int const n = run(&sy,1u<<TICK_WWV | 1u<<TICK_WWVH,3,-1,NULL);
  assert(n == 3);
  for(int i=0; i < n; i++){
    struct tick_second const *r = &Rep[i];
    assert(r->complete);
    int64_t const sec = S0 + (int64_t)i * SR;
    assert(r->second_index == sec);
    assert(r->window_start_index == sec - 240);
    assert(r->step == 2 && r->nprofile == 840);
    double const wwv_true = sec + 4.2e-3 * SR, wwvh_true = sec + 29.0e-3 * SR;
    double const tol = 0.02e-3 * SR;
    printf("  second %d: wwv %+.4f wwvh %+.4f samples from truth; snr %.1f %.1f dB\n",i,
           r->peak_index[TICK_WWV] - wwv_true,r->peak_index[TICK_WWVH] - wwvh_true,
           r->peak_snr_db[TICK_WWV],r->peak_snr_db[TICK_WWVH]);
    assert(fabs(r->peak_index[TICK_WWV] - wwv_true) <= tol);
    assert(fabs(r->peak_index[TICK_WWVH] - wwvh_true) <= tol);
    // At WWV's lag the WWVH template sees WWV's burst aligned: orthogonal
    int const k = (int)lround((wwv_true - r->window_start_index) / r->step);
    double const leak = cabsf(r->profile[TICK_WWVH][k]) / maxmag(r->profile[TICK_WWVH],r->nprofile);
    printf("  WWVH template at WWV's lag: %.1f dB\n",20*log10(leak));
    assert(20*log10(leak) <= -20);
    // and the other way round
    int const kh = (int)lround((wwvh_true - r->window_start_index) / r->step);
    double const leak2 = cabsf(r->profile[TICK_WWV][kh]) / maxmag(r->profile[TICK_WWV],r->nprofile);
    assert(20*log10(leak2) <= -20);
  }
  printf("test_two_stations_separate passed\n");
}

// The absent station keeps its profile and reports a low SNR. The brief's
// "< 6 dB" cannot hold under its own SNR (peak power over median power):
// noise alone scores about 10 dB, and a misaligned 1000 Hz burst leaks into
// the 1200 Hz template about 10 dB below WWV's own peak. So "low" here means
// at least 6 dB under the present station's
static void test_absent_station_reported(void){
  struct synth sy = { .tick_ms = {4.2,0}, .tick_amp = {1.0,0}, .noise_rms = 0.3,
                      .first_som = 10, .phase = -1.1 };
  int const n = run(&sy,1u<<TICK_WWV | 1u<<TICK_WWVH,3,-1,NULL);
  assert(n == 3);
  for(int i=0; i < n; i++){
    struct tick_second const *r = &Rep[i];
    assert(r->complete);
    int nonzero = 0;
    for(int k=0; k < r->nprofile; k++){
      assert(isfinite(crealf(r->profile[TICK_WWVH][k])) && isfinite(cimagf(r->profile[TICK_WWVH][k])));
      nonzero += r->profile[TICK_WWVH][k] != 0;
    }
    assert(nonzero == r->nprofile);
    assert(isfinite(r->peak_index[TICK_WWVH]) && isfinite(r->peak_snr_db[TICK_WWVH]));
    printf("  second %d: snr wwv %.1f dB, absent wwvh %.1f dB\n",i,
           r->peak_snr_db[TICK_WWV],r->peak_snr_db[TICK_WWVH]);
    assert(r->peak_snr_db[TICK_WWV] > 20);
    assert(r->peak_snr_db[TICK_WWVH] < r->peak_snr_db[TICK_WWV] - 6);
  }
  printf("test_absent_station_reported passed\n");
}

static void test_skipped_second_still_reports(void){
  // Whole seconds 28, 29, 30: 29 carries no tick
  struct synth sy = { .tick_ms = {4.2,29.0}, .tick_amp = {1.0,0.5}, .noise_rms = 0.05,
                      .first_som = 28, .phase = 2.0 };
  unsigned const all = 1u<<TICK_WWV | 1u<<TICK_WWVH | 1u<<TICK_BPM;
  int const n = run(&sy,all,3,-1,NULL);
  assert(n == 3);
  struct tick_second const *r = &Rep[1];   // second 29
  assert(r->complete);
  for(int t=0; t < TICK_NTEMPLATES; t++){
    double const m = maxmag(r->profile[t],r->nprofile);
    assert(m > 0 && isfinite(m));
    assert(isfinite(r->peak_index[t]) && isfinite(r->peak_snr_db[t]));
    assert(r->peak_index[t] >= r->window_start_index
           && r->peak_index[t] <= r->window_start_index + (double)(r->nprofile-1) * r->step);
    // and the tick is gone: second 29's profile peak sits far below second 28's
    assert(m < 0.1 * maxmag(Rep[0].profile[t],Rep[0].nprofile));
  }
  printf("  second 29: snr wwv %.1f wwvh %.1f bpm %.1f dB\n",
         r->peak_snr_db[TICK_WWV],r->peak_snr_db[TICK_WWVH],r->peak_snr_db[TICK_BPM]);
  printf("test_skipped_second_still_reports passed\n");
}

static void test_envelope_rate(void){
  struct synth sy = { .tone_hz = 100, .tone_depth = 0.5, .first_som = 10, .phase = 0.3 };
  int const n = run(&sy,1u<<TICK_WWV,2,-1,NULL);
  assert(n == 2);
  for(int i=0; i < n; i++){
    struct tick_second const *r = &Rep[i];
    assert(r->complete && r->envelope_scale > 0);
    double mean = 0;
    double complex bin = 0, off = 0;
    for(int k=0; k < TICK_ENVELOPE; k++){
      double const e = r->envelope[k] * (double)r->envelope_scale;
      mean += e;
      bin += e * cexp(-I*2*M_PI*100*k/1000.0);
      off += e * cexp(-I*2*M_PI*150*k/1000.0);
    }
    mean /= TICK_ENVELOPE;
    double const amp = 2 * cabs(bin) / TICK_ENVELOPE;
    double const amp150 = 2 * cabs(off) / TICK_ENVELOPE;
    printf("  second %d: mean %.5f, 100 Hz amplitude %.5f (150 Hz %.2g)\n",i,mean,amp,amp150);
    assert(fabs(mean - 1) < 0.02);
    assert(fabs(amp - 0.5) <= 0.02 * 0.5);
    assert(amp150 < 0.01);
    // The 100 Hz phase follows the sample clock: envelope[k] sits at second + 24 k
    int64_t const sec = S0 + (int64_t)i * SR;
    double const want = fmod(2*M_PI*100*(double)sec/SR,2*M_PI);
    double const got = carg(bin);
    double d = remainder(got - want,2*M_PI);
    assert(fabs(d) < 0.05);
  }
  printf("test_envelope_rate passed\n");
}

static void test_marker_energy(void){
  struct synth sy = { .tick_ms = {4.2,29.0}, .tick_amp = {1.0,1.0}, .marker_hz = 1500, .marker_amp = 0.5,
                      .noise_rms = 0.02, .first_som = 0, .phase = 1.3 };
  int const n = run(&sy,1u<<TICK_WWV | 1u<<TICK_WWVH,2,-1,NULL);
  assert(n == 2);
  struct tick_second const *r = &Rep[0];
  assert(r->complete);
  printf("  second 0 markers: 1000 %.3g, 1200 %.3g, 1500 %.3g\n",
         r->marker_energy[0],r->marker_energy[1],r->marker_energy[2]);
  // AM depth 0.5: each sideband 0.25, energy = their summed power 0.125
  assert(fabs(r->marker_energy[2] - 0.125) < 0.125 * 0.01);
  assert(10*log10(r->marker_energy[2] / r->marker_energy[0]) >= 20);
  assert(10*log10(r->marker_energy[2] / r->marker_energy[1]) >= 20);
  // second 1 has no marker
  assert(Rep[1].marker_energy[2] < 1e-3 * r->marker_energy[2]);
  printf("test_marker_energy passed\n");
}

// A missing block returns -1 once; the report whose span held it says
// incomplete, and the next second reports complete again
static void test_gap_marks_incomplete(void){
  struct synth sy = { .tick_ms = {4.2,0}, .tick_amp = {1.0,0}, .noise_rms = 0.02, .first_som = 10 };
  int gaps;
  int const skip = (S0 + SR + SR/2) / BATCH;   // mid second 1
  int const n = run(&sy,1u<<TICK_WWV,3,skip,&gaps);
  assert(gaps == 1);
  assert(n == 3);
  assert(Rep[0].complete);
  assert(!Rep[1].complete);
  assert(Rep[2].complete);
  printf("test_gap_marks_incomplete passed\n");
}

int main(void){
  test_two_stations_separate();
  test_absent_station_reported();
  test_skipped_second_still_reports();
  test_envelope_rate();
  test_marker_energy();
  test_gap_marks_incomplete();
  printf("all tick tests passed\n");
  exit(0);
}

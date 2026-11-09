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
#include <time.h>
#include "tick.h"
#include "tick_feed.h"
#include "pilot.h"
#include "tick_status.h"

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
    for(int k=0; k < TICK_ENVELOPE_LEN; k++){
      double const e = r->envelope[k] * (double)r->envelope_scale;
      mean += e;
      bin += e * cexp(-I*2*M_PI*100*k/1000.0);
      off += e * cexp(-I*2*M_PI*150*k/1000.0);
    }
    mean /= TICK_ENVELOPE_LEN;
    double const amp = 2 * cabs(bin) / TICK_ENVELOPE_LEN;
    double const amp150 = 2 * cabs(off) / TICK_ENVELOPE_LEN;
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

// ---- radiod glue (tick_feed.c): naming the second, placing the window, epochs ----

static void test_second_index_math(void){
  int64_t const S = 1400000000LL * 1000000000LL;
  // 24 kHz, reference index 1 000 000 = UTC S + 0.25 s: second S+1 sits 0.75 s later
  int64_t const k = tick_second_index(S + 1000000000LL,1000000,S + 250000000LL,24000);
  printf("  S+1 at index %lld\n",(long long)k);
  assert(k == 1018000);
  assert(tick_second_index(S,1000000,S + 250000000LL,24000) == 1000000 - 6000);
  // A reference 100 s back, at an input rate, stays exact to the sample
  assert(tick_second_index(S + 100000000000LL,7,S,64800000.0) == 7 + 6480000000LL);
  printf("test_second_index_math passed\n");
}

static void test_parse_templates(void){
  assert(tick_parse_templates("wwv") == 1u << TICK_WWV);
  assert(tick_parse_templates("wwv+wwvh") == ((1u << TICK_WWV) | (1u << TICK_WWVH)));
  assert(tick_parse_templates("WWV+wwvh+bpm") == 7u);
  assert(tick_parse_templates("wwv+bpm") == ((1u << TICK_WWV) | (1u << TICK_BPM)));
  assert(tick_parse_templates("") == 0);
  assert(tick_parse_templates("wwv+chu") == 0);
  assert(tick_parse_templates("wwv++wwvh") == 0);
  assert(tick_parse_templates(NULL) == 0);
  printf("test_parse_templates passed\n");
}

#define DECIM 4
#define IN_SR (SR * DECIM)
#define BASE_IN 4000000004LL           // input index of the first block
#define UTC0 (1400000000LL * 1000000000LL) // UTC of the second at output index S0

struct feed_opts {
  int64_t ref_shift_ns;   // reference names UTC this much later than truth
  int epoch_block;        // block where epoch 2 begins (-1: none)
  bool straddle;          // that block straddles the loss instead
  int pilot_from_block;   // from this block on, the reference counts as the pilot's (-1: never)
  int64_t pilot_shift_ns; // the pilot reference's own error
  double delay_out;       // tick channel filter delay, output samples
  int64_t loss_ns;        // UTC lost at the epoch change: the input counter stays
                          // contiguous, so from epoch 2 on the reference names UTC this much later
};
static struct tick_meta Meta[MAXREP];

static int run_feed(struct synth const *sy,unsigned mask,int seconds,struct feed_opts const *o){
  struct tick_config const cfg = { .samprate = SR, .templates = mask };
  struct tick_state *t = tick_create(&cfg);
  assert(t != NULL);
  struct tick_feed f;
  tick_feed_init(&f,DECIM,IN_SR,SR,o->delay_out);
  static float complex buf[BATCH];
  int nrep = 0;
  int64_t const total = S0 + (int64_t)(seconds-1) * SR + tick_span_end(t) + BATCH;
  int block = 0;
  for(int64_t index = 0; index < total; index += BATCH, block++){
    for(int i=0; i < BATCH; i++)
      buf[i] = synth_sample(sy,index + i);
    uint64_t const epoch = (o->epoch_block >= 0 && block >= o->epoch_block + (o->straddle ? 1 : 0)) ? 2 : 1;
    bool const straddle = o->straddle && block == o->epoch_block;
    bool const pilot = o->pilot_from_block >= 0 && block >= o->pilot_from_block;
    // The UTC second at output index S0 enters at input BASE_IN + DECIM*S0
    struct tick_ref const ref = {
      .input_pos = (double)(BASE_IN + DECIM * (int64_t)S0),
      .utc_ns = UTC0 + (pilot ? o->pilot_shift_ns : o->ref_shift_ns) + (epoch == 2 ? o->loss_ns : 0),
      .epoch = epoch,
    };
    struct tick_second r;
    struct tick_meta m;
    int const ret = tick_feed_block(t,&f,buf,BATCH,BASE_IN + DECIM * index,epoch,straddle,&ref,pilot,&r,&m);
    if(ret == 1 && nrep < MAXREP){
      Rep[nrep] = r;
      Meta[nrep] = m;
      nrep++;
    }
  }
  tick_destroy(t);
  return nrep;
}

static void test_coarse_nominal_second(void){
  struct synth sy = { .tick_ms = {3.0,0}, .tick_amp = {0.5,0}, .noise_rms = 0.01, .first_som = 10, .phase = 0.7 };
  struct feed_opts o = { .epoch_block = -1, .pilot_from_block = -1 };
  Rng = 1;
  int const n0 = run_feed(&sy,1u << TICK_WWV,3,&o);
  struct tick_second a[3];
  struct tick_meta am[3];
  assert(n0 >= 3);
  memcpy(a,Rep,sizeof(a));
  memcpy(am,Meta,sizeof(am));
  o.ref_shift_ns = 5000000; // the reference says UTC is 5 ms later than it is
  Rng = 1;
  int const n1 = run_feed(&sy,1u << TICK_WWV,3,&o);
  assert(n1 == n0);
  for(int i=0; i < 3; i++){
    struct tick_second const *b = &Rep[i];
    double const truth = S0 + (double)i * SR + 3.0e-3 * SR;
    printf("  second %d: exact ref peak %.4f, 5 ms ref peak %.4f (truth %.1f); window starts %lld vs %lld\n",
           i,a[i].peak_index[TICK_WWV],b->peak_index[TICK_WWV],truth,
           (long long)a[i].window_start_index,(long long)b->window_start_index);
    assert(a[i].complete && b->complete);
    assert(a[i].second_index == S0 + (int64_t)i * SR);
    assert(b->second_index == a[i].second_index - 120); // 5 ms at 24 kHz
    assert(am[i].second_utc_ns == UTC0 + (int64_t)i * 1000000000LL);
    assert(Meta[i].second_utc_ns == am[i].second_utc_ns); // same second named
    assert(!Meta[i].used_pilot && Meta[i].epoch == 1);
    assert(fabs(b->peak_index[TICK_WWV] - a[i].peak_index[TICK_WWV]) < 0.01);
    assert(fabs(a[i].peak_index[TICK_WWV] - truth) < 0.5);
  }
  // The tick channel's own filter delay moves the window, never the measured peak
  o.ref_shift_ns = 0;
  o.delay_out = 60;
  Rng = 1;
  assert(run_feed(&sy,1u << TICK_WWV,3,&o) >= 3);
  assert(Rep[1].second_index == a[1].second_index + 60);
  assert(fabs(Rep[1].peak_index[TICK_WWV] - a[1].peak_index[TICK_WWV]) < 0.01);
  printf("test_coarse_nominal_second passed\n");
}

static void test_window_spanning_epoch_is_flagged(void){
  struct synth sy = { .tick_ms = {3.0,0}, .tick_amp = {0.5,0}, .noise_rms = 0.01, .first_som = 10 };
  // Second 2 starts at output S0 + 2 SR; its window runs -10..+60 ms. Block 125
  // starts at 125*480 = 60000 = S0 + 2 SR + 0 ms; block 126 at +20 ms: inside the window
  int const blk = (int)((S0 + 2 * SR) / BATCH) + 1;
  for(int straddle = 0; straddle <= 1; straddle++){
    struct feed_opts const o = { .epoch_block = blk, .straddle = straddle, .pilot_from_block = -1 };
    Rng = 2;
    int const n = run_feed(&sy,1u << TICK_WWV,5,&o);
    assert(n >= 5);
    for(int i=0; i < 5; i++)
      printf("  %s: second %d complete %d epoch %llu peak %.3f\n",straddle ? "straddling block" : "block boundary",
             i,Rep[i].complete,(unsigned long long)Meta[i].epoch,Rep[i].peak_index[TICK_WWV]);
    assert(Rep[1].complete && Meta[1].epoch == 1);
    assert(Rep[2].second_index == S0 + 2 * SR);
    assert(!Rep[2].complete && isnan(Rep[2].peak_index[TICK_WWV]));
    // Second 3 was named (at S0 + 2 SR - 10 ms) in epoch 1 and completes in epoch 2
    assert(Rep[3].second_index == S0 + 3 * SR);
    assert(!Rep[3].complete && isnan(Rep[3].peak_index[TICK_WWV]) && Meta[3].epoch == 2);
    assert(Rep[4].complete && Meta[4].epoch == 2);
  }
  printf("test_window_spanning_epoch_is_flagged passed\n");
}

// A loss keeps the input counter contiguous while UTC jumps by the lost time:
// a second named on the old mapping must not report complete on new data,
// and the next second follows the new mapping
static void test_loss_moves_the_mapping(void){
  struct synth sy = { .tick_ms = {3.0,0}, .tick_amp = {0.5,0}, .noise_rms = 0.01, .first_som = 10 };
  int const blk = (int)((S0 + 2 * SR) / BATCH) + 1;
  int64_t const L = 300000000; // 300 ms lost
  struct feed_opts const o = { .epoch_block = blk, .pilot_from_block = -1, .loss_ns = L };
  Rng = 4;
  int const n = run_feed(&sy,1u << TICK_WWV,6,&o);
  assert(n >= 5);
  for(int i=0; i < 5; i++)
    printf("  second %d: index %lld complete %d epoch %llu utc +%lld ms\n",i,(long long)Rep[i].second_index,
           Rep[i].complete,(unsigned long long)Meta[i].epoch,(long long)((Meta[i].second_utc_ns - UTC0) / 1000000));
  int64_t const shift = L / 1000000 * SR / 1000; // 7200 samples
  assert(Rep[3].second_index == S0 + 3 * SR);   // named on the old mapping
  assert(!Rep[3].complete && isnan(Rep[3].peak_index[TICK_WWV]));
  assert(Rep[4].second_index == S0 + 4 * SR - shift);
  assert(Meta[4].second_utc_ns == UTC0 + 4000000000LL && Meta[4].epoch == 2);
  assert(Rep[4].complete);
  printf("test_loss_moves_the_mapping passed\n");
}

// The pilot's raw edge runs late by the pilot channel's filter delay: the
// reference takes it off
static void test_ref_from_pilot(void){
  struct pilot_report r = { .valid = true, .locked = true, .epoch = 9,
                            .edge_input_sample = 1000000000LL, .edge_input_frac = 0.25,
                            .edge_utc_ns = UTC0, .filter_delay_s = 0.0025 };
  double const in = 64.8e6;
  struct tick_ref ref;
  assert(tick_ref_from_pilot(&r,in,&ref));
  printf("  edge %lld + 0.25, delay 2.5 ms -> input_pos %.2f\n",(long long)r.edge_input_sample,ref.input_pos);
  assert(fabs(ref.input_pos - (1000000000.25 - 162000)) < 1e-3);
  assert(ref.utc_ns == UTC0 && ref.epoch == 9);
  r.locked = false;
  assert(!tick_ref_from_pilot(&r,in,&ref));
  r.locked = true;
  r.valid = false;
  assert(!tick_ref_from_pilot(&r,in,&ref));
  printf("test_ref_from_pilot passed\n");
}

// A pilot reference arriving mid-run: the second already named on the host
// clock keeps that alignment (tick.c's half-second rule) and says so; the
// next second follows the pilot
static void test_used_pilot_follows_the_named_second(void){
  struct synth sy = { .tick_ms = {3.0,0}, .tick_amp = {0.5,0}, .noise_rms = 0.01, .first_som = 10 };
  // Host clock 3 ms off; pilot exact, from block 60 (output 28800 = S0 + 0.7 s)
  struct feed_opts const o = { .epoch_block = -1, .ref_shift_ns = 3000000, .pilot_from_block = 60, .pilot_shift_ns = 0 };
  Rng = 3;
  int const n = run_feed(&sy,1u << TICK_WWV,4,&o);
  assert(n >= 4);
  for(int i=0; i < 4; i++)
    printf("  second %d: index %lld used_pilot %d utc +%lld ms\n",i,(long long)Rep[i].second_index,
           Meta[i].used_pilot,(long long)((Meta[i].second_utc_ns - UTC0) / 1000000));
  assert(!Meta[0].used_pilot && Rep[0].second_index == S0 - 72);
  // Second 1 was named ~1 s ahead (S0 - 10 ms), on the host clock
  assert(!Meta[1].used_pilot && Rep[1].second_index == S0 + SR - 72);
  assert(Meta[2].used_pilot && Rep[2].second_index == S0 + 2 * SR);
  assert(Meta[3].used_pilot && Meta[3].second_utc_ns == UTC0 + 3000000000LL);
  printf("test_used_pilot_follows_the_named_second passed\n");
}

static void test_pilot_slot(void){
  static struct tick_pilot_slot slot = TICK_PILOT_SLOT_INIT;
  struct tick_ref r;
  assert(!tick_slot_read(&slot,&r)); // nothing published yet
  struct tick_ref const p = { .input_pos = 123456789012.0, .utc_ns = UTC0, .epoch = 7 };
  tick_slot_publish(&slot,&p);
  assert(tick_slot_read(&slot,&r));
  assert(r.input_pos == p.input_pos && r.utc_ns == p.utc_ns && r.epoch == 7);
  double const in = 64.8e6;
  assert(tick_ref_fresh(&r,7,(int64_t)(p.input_pos + 119 * in),in));
  assert(!tick_ref_fresh(&r,7,(int64_t)(p.input_pos + 121 * in),in)); // older than 120 s
  assert(!tick_ref_fresh(&r,8,(int64_t)p.input_pos + 1000,in));      // another epoch
  printf("test_pilot_slot passed\n");
}

// Correlation cost on this machine: three templates, markers and envelope, per second
static void test_report_cost(void){
  struct tick_config const cfg = { .samprate = SR, .templates = 7 };
  struct tick_state *t = tick_create(&cfg);
  assert(t != NULL);
  struct synth sy = { .tick_ms = {3.0,0}, .tick_amp = {0.5,0}, .noise_rms = 0.1, .first_som = 10 };
  int const seconds = 30;
  int64_t const total = (int64_t)seconds * SR;
  float complex *all = malloc(total * sizeof(*all));
  assert(all != NULL);
  for(int64_t i=0; i < total; i++)
    all[i] = synth_sample(&sy,i);
  struct timespec a,b;
  int nrep = 0;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&a);
  for(int64_t index = 0; index + BATCH <= total; index += BATCH){
    struct tick_second r;
    if(tick_process(t,all + index,BATCH,index,second_of(index),&r) == 1)
      nrep++;
  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&b);
  double const ms = ((b.tv_sec - a.tv_sec) * 1e9 + (b.tv_nsec - a.tv_nsec)) / 1e6;
  printf("  %d reports, %.2f ms CPU per second of 24 kHz input (all three templates)\n",nrep,ms / seconds);
  assert(nrep >= seconds - 2);
  free(all);
  tick_destroy(t);
  printf("test_report_cost passed\n");
}

// A filled report through encode_tick_status and decode_tick_status: every
// field comes back, and the packet stays inside what one datagram holds
// Accepted size of one tick group with three 840-value profiles and the envelope
// (measured 12,176 bytes; three int16-pair profiles are 10 KB alone). Above the
// 9 KB first hoped for; the group is sent last in the packet for that reason
#define TICK_GROUP_BYTES_ACCEPTED 12500

static void decode_group(struct tick_status *d,uint8_t const *cp,uint8_t const *end){
  while(cp < end){
    int const type = *cp++;
    int optlen = *cp++;
    if(optlen & 0x80){
      int const n = optlen & 0x7f;
      optlen = 0;
      for(int i=0; i < n; i++)
        optlen = (optlen << 8) | *cp++;
    }
    assert(decode_tick_status(d,(enum status_type)type,cp,optlen));
    cp += optlen;
  }
}

static void test_tick_once_only(void){
  uint64_t sent[2] = {0,0};
  assert(!tick_status_take(0,&sent[0]));          // no report yet
  assert(tick_status_take(1,&sent[0]));           // first packet carries report 1
  assert(!tick_status_take(1,&sent[0]));          // later packets do not
  assert(!tick_status_take(1,&sent[0]));
  assert(tick_status_take(1,&sent[1]));           // the other destination still gets it once
  assert(!tick_status_take(1,&sent[1]));
  assert(tick_status_take(2,&sent[0]));           // a new report goes out again
  assert(tick_status_take(3,&sent[0]));           // a skipped number is fine
  printf("test_tick_once_only passed\n");
}

static void test_tick_tlv_roundtrip(void){
  struct tick_second *s = calloc(1,sizeof(*s));
  struct tick_status *d = calloc(1,sizeof(*d));
  assert(s != NULL && d != NULL);
  int const decim = 2560;      // 64.8 MHz / 24 kHz is 2700; any integer exercises the arithmetic
  int64_t const base = 123456789012LL;
  s->second_index = 5000000; s->window_start_index = 5000000 - 240; s->step = 2; s->nprofile = 840;
  for(int k=0; k < TICK_NTEMPLATES; k++){
    for(int i=0; i < s->nprofile; i++)
      s->profile[k][i] = (float)(37.5 * (k+1) * sin(0.01*i)) + I * (float)(-12.25 * cos(0.02*i*(k+1)));
    s->peak_index[k] = s->window_start_index + 100.25 + 7 * k;
    s->peak_snr_db[k] = 20.5f + k;
  }
  for(int m=0; m < TICK_MARKER_TONES; m++)
    s->marker_energy[m] = 0.125 * (m+1);
  for(int i=0; i < TICK_ENVELOPE_LEN; i++)
    s->envelope[i] = (int16_t)(i * 31 - 15000);
  s->envelope_scale = 0.00042f;
  s->complete = true;
  uint8_t pkt[16384];
  uint8_t *bp = pkt;
  encode_tick_status(&bp,s,1700000000000000000LL,true,7,base,decim);
  size_t const len = bp - pkt;
  printf("  tick group %zu bytes (three profiles of %d, envelope)\n",len,s->nprofile);
  assert(len < TICK_GROUP_BYTES_ACCEPTED);

  decode_group(d,pkt,bp);
  assert(d->valid && d->second_utc_ns == 1700000000000000000LL);
  assert(d->window_start_input == base + (int64_t)s->window_start_index * decim);
  assert(d->step_input == 2 * decim);
  assert(d->complete && d->used_pilot && d->templates == 7);
  for(int k=0; k < TICK_NTEMPLATES; k++){
    assert(d->nprofile[k] == 840);
    float want = 0;   // the transmitted scale is the largest |re| or |im| of the profile
    for(int i=0; i < 840; i++){
      want = fmaxf(want,fabsf(crealf(s->profile[k][i])));
      want = fmaxf(want,fabsf(cimagf(s->profile[k][i])));
    }
    assert(d->profile_scale[k] == want);
    float const tol = d->profile_scale[k] / 32767.0f * 1.01f;
    for(int i=0; i < 840; i++){
      assert(fabsf(crealf(d->profile[k][i]) - crealf(s->profile[k][i])) <= tol);
      assert(fabsf(cimagf(d->profile[k][i]) - cimagf(s->profile[k][i])) <= tol);
    }
    assert(d->have_peaks);
    assert(d->peak_offset[k] == (float)((100.25 + 7 * k) * decim));
    assert(d->peak_snr_db[k] == (float)s->peak_snr_db[k]);
  }
  assert(d->have_markers);
  for(int m=0; m < TICK_MARKER_TONES; m++)
    assert(d->marker_energy[m] == (float)s->marker_energy[m]);
  assert(d->nenvelope == TICK_ENVELOPE_LEN && d->envelope_scale == s->envelope_scale);
  for(int i=0; i < TICK_ENVELOPE_LEN; i++)
    assert(d->envelope[i] == s->envelope[i]);

  // An incomplete report: zeros, NaN peaks, bit0 clear, still delivered
  tick_second_incomplete(s);
  s->complete = false;
  bp = pkt;
  encode_tick_status(&bp,s,1700000001000000000LL,false,1u << TICK_WWV,base,decim);
  memset(d,0,sizeof(*d));
  decode_group(d,pkt,bp);
  assert(d->valid && !d->complete && !d->used_pilot && d->templates == 1);
  assert(d->nprofile[TICK_WWV] == 840 && d->nprofile[TICK_WWVH] == 0);
  assert(maxmag(d->profile[TICK_WWV],840) == 0);
  assert(d->have_peaks && isnan(d->peak_offset[TICK_WWV]) && isnan(d->peak_offset[TICK_WWVH]));
  // A non-contiguous mask (wwv + bpm): wwvh's slot stays empty, its peak NaN
  s->complete = true;
  for(int i=0; i < 840; i++)
    s->profile[TICK_WWV][i] = s->profile[TICK_BPM][i] = (float)i - 400.0f + I * 3.0f;
  s->profile[TICK_WWV][5] = NAN;            // a non-finite value must not poison the scale
  for(int k=0; k < TICK_NTEMPLATES; k++){ s->peak_index[k] = s->window_start_index + 10; s->peak_snr_db[k] = 9; }
  bp = pkt;
  encode_tick_status(&bp,s,1700000002000000000LL,false,5,base,decim);
  memset(d,0,sizeof(*d));
  decode_group(d,pkt,bp);
  assert(d->templates == 5 && d->nprofile[TICK_WWV] == 840 && d->nprofile[TICK_BPM] == 840);
  assert(d->nprofile[TICK_WWVH] == 0);
  assert(isfinite(d->profile_scale[TICK_WWV]) && d->profile_scale[TICK_WWV] == 439.0f);
  assert(crealf(d->profile[TICK_WWV][5]) == 0);
  assert(isnan(d->peak_offset[TICK_WWVH]) && d->peak_offset[TICK_BPM] == 10.0f * decim);
  free(s); free(d);
  printf("test_tick_tlv_roundtrip passed\n");
}

int main(void){
  test_two_stations_separate();
  test_absent_station_reported();
  test_skipped_second_still_reports();
  test_envelope_rate();
  test_marker_energy();
  test_gap_marks_incomplete();
  test_second_index_math();
  test_parse_templates();
  test_coarse_nominal_second();
  test_window_spanning_epoch_is_flagged();
  test_used_pilot_follows_the_named_second();
  test_loss_moves_the_mapping();
  test_ref_from_pilot();
  test_pilot_slot();
  test_report_cost();
  test_tick_tlv_roundtrip();
  test_tick_once_only();
  printf("all tick tests passed\n");
  exit(0);
}

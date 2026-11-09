// Time-station tick detector: what radiod measures each second on a WWV/WWVH/BPM channel
//
// For each enabled template, a complex profile: the carrier-centred IQ
// correlated against a real template cos(2 pi f k/fs) over the burst,
//   profile[i] = sum_k s[w + i*step + k] cos(2 pi f k/fs),  k = 0..L-1,
// for lags w + i*step across -10..+60 ms of the nominal second. AM puts a tick
// on both sidebands at +/-f; the real template collects both, and the
// profile's phase follows the carrier. Then the peak of |profile| with a
// parabola through it and its neighbours, and an SNR: peak power over the
// median power of the lags more than 3 ms from the peak.
// Per marker tone (1000, 1200, 1500 Hz), the energy over 0..800 ms by
// Goertzel; one recurrence serves the bins at +f and -f, and the energy is
// (|X(+f)|^2 + |X(-f)|^2)/N^2, the mean square of the modulating tone.
// The AM envelope: |s| through a Hamming-windowed sinc low-pass (250 Hz),
// one output per fs/1000 samples, envelope[i] centred on second + i*fs/1000.
//
// Samples go into a ring as they arrive. Each nominal second the caller
// names waits in a short queue until the stream passes the end of its span
// (about second + 1.008 s, the envelope's last filter tap); then its report
// is computed from a linear copy of the span. A report whose span does not
// lie wholly in contiguous ring data comes out with complete = false and
// zeros (peak_index NaN). Nothing here judges which station is present.
// Copyright 2026, HamSCI
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "tick.h"

#define WINDOW_BEGIN_S (-0.010)
#define WINDOW_LEN_S 0.070
#define STEP_S 0.0001
#define SNR_EXCLUDE_S 0.003
#define MARKER_LEN_S 0.8
#define ENV_RATE 1000
#define ENV_CUTOFF_HZ 250.0
#define ENV_HALF_OUT 8        // filter half-length in envelope samples
#define PENDING 4
#define RING_SECONDS 2

static double const Tmpl_hz[TICK_NTEMPLATES] = {1000, 1200, 1000};
static double const Tmpl_cycles[TICK_NTEMPLATES] = {5, 6, 10};
static char const * const Tmpl_name[TICK_NTEMPLATES] = {"wwv", "wwvh", "bpm"};
static double const Marker_hz[TICK_MARKERS] = {1000, 1200, 1500};

struct tick_state {
  struct tick_config cfg;
  int step, nprofile, win_off;   // win_off: first lag relative to the nominal second
  int excl;                      // SNR exclusion half-width, samples
  int decim, env_half, env_ntaps;
  int marker_n;
  int tlen[TICK_NTEMPLATES];
  float *tmpl[TICK_NTEMPLATES];
  double *env_taps;
  int span_begin, span_end;      // report span [second + span_begin, second + span_end)
  int span_len;

  float complex *ring;
  int ringsize;
  bool have_end;
  int64_t end;                   // counter of the next sample expected
  int64_t valid_from;            // ring data contiguous from here to end

  int64_t pending[PENDING];
  int npending;
  bool have_announced;
  int64_t last_announced;

  // Scratch, allocated once
  float complex *lin;            // the span, linear
  float *mag;                    // |lin|
  double *power;                 // median work
};

char const *tick_template_name(enum tick_template k){
  return (k >= 0 && k < TICK_NTEMPLATES) ? Tmpl_name[k] : "?";
}

int tick_span_end(struct tick_state const *t){
  return t->span_end;
}

struct tick_state *tick_create(struct tick_config const *cfg){
  if(cfg == NULL || cfg->samprate < 4000 || cfg->samprate % ENV_RATE != 0
     || cfg->templates == 0 || (cfg->templates >> TICK_NTEMPLATES) != 0)
    return NULL;
  struct tick_state *t = calloc(1,sizeof(*t));
  if(t == NULL)
    return NULL;
  t->cfg = *cfg;
  double const fs = cfg->samprate;
  t->step = (int)lround(STEP_S * fs);
  if(t->step < 1)
    t->step = 1;
  t->win_off = -(int)lround(-WINDOW_BEGIN_S * fs);
  t->nprofile = (int)lround(WINDOW_LEN_S * fs) / t->step;
  t->excl = (int)lround(SNR_EXCLUDE_S * fs);
  t->decim = cfg->samprate / ENV_RATE;
  t->env_half = ENV_HALF_OUT * t->decim;
  t->env_ntaps = 2 * t->env_half + 1;
  t->marker_n = (int)lround(MARKER_LEN_S * fs);
  if(t->nprofile > TICK_MAX_PROFILE){
    free(t);
    return NULL;
  }
  bool ok = true;
  int maxlen = 0;
  for(int k=0; k < TICK_NTEMPLATES; k++){
    t->tlen[k] = (int)lround(Tmpl_cycles[k] * fs / Tmpl_hz[k]);
    if(t->tlen[k] > maxlen)
      maxlen = t->tlen[k];
    t->tmpl[k] = malloc(t->tlen[k] * sizeof(float));
    ok = ok && t->tmpl[k];
    for(int j=0; ok && j < t->tlen[k]; j++)
      t->tmpl[k][j] = (float)cos(2*M_PI*Tmpl_hz[k]*j/fs);
  }
  // Hamming-windowed sinc, unity gain at DC
  t->env_taps = malloc(t->env_ntaps * sizeof(double));
  if(t->env_taps != NULL){
    double sum = 0;
    for(int k=0; k < t->env_ntaps; k++){
      double const x = 2 * ENV_CUTOFF_HZ / fs * (k - t->env_half);
      double const sinc = x == 0 ? 1.0 : sin(M_PI*x) / (M_PI*x);
      double const w = 0.54 - 0.46 * cos(2*M_PI*k/(t->env_ntaps - 1));
      t->env_taps[k] = sinc * w;
      sum += t->env_taps[k];
    }
    for(int k=0; k < t->env_ntaps; k++)
      t->env_taps[k] /= sum;
  }
  // Span: the window's lags plus the longest template, the marker interval,
  // and the envelope's filter taps
  t->span_begin = t->win_off < -t->env_half ? t->win_off : -t->env_half;
  int e = t->win_off + (t->nprofile - 1) * t->step + maxlen;
  if(t->marker_n > e)
    e = t->marker_n;
  int const env_last = (TICK_ENVELOPE - 1) * t->decim + t->env_half + 1;
  if(env_last > e)
    e = env_last;
  t->span_end = e;
  t->span_len = t->span_end - t->span_begin;
  t->ringsize = RING_SECONDS * cfg->samprate;
  if(t->span_len > t->ringsize)
    ok = false;
  t->ring = calloc(t->ringsize,sizeof(*t->ring));
  t->lin = malloc(t->span_len * sizeof(*t->lin));
  t->mag = malloc(t->span_len * sizeof(*t->mag));
  t->power = malloc(TICK_MAX_PROFILE * sizeof(*t->power));
  if(!ok || !t->env_taps || !t->ring || !t->lin || !t->mag || !t->power){
    tick_destroy(t);
    return NULL;
  }
  tick_reset(t);
  return t;
}

void tick_destroy(struct tick_state *t){
  if(t == NULL)
    return;
  for(int k=0; k < TICK_NTEMPLATES; k++)
    free(t->tmpl[k]);
  free(t->env_taps);
  free(t->ring);
  free(t->lin);
  free(t->mag);
  free(t->power);
  free(t);
}

void tick_reset(struct tick_state *t){
  if(t == NULL)
    return;
  t->have_end = false;
  t->end = 0;
  t->valid_from = 0;
  t->npending = 0;
  t->have_announced = false;
  t->last_announced = 0;
}

void tick_break(struct tick_state *t,int64_t index){
  if(t != NULL && t->have_end && index > t->valid_from)
    t->valid_from = index;
}

static int ring_pos(struct tick_state const *t,int64_t index){
  int64_t r = index % t->ringsize;
  if(r < 0)
    r += t->ringsize;
  return (int)r;
}

// k-th smallest of a[0..n-1] (0-based); reorders a
static double select_kth(double *a,int n,int k){
  int lo = 0, hi = n - 1;
  while(lo < hi){
    double const pivot = a[(lo + hi) / 2];
    int i = lo, j = hi;
    while(i <= j){
      while(a[i] < pivot)
        i++;
      while(a[j] > pivot)
        j--;
      if(i <= j){
        double const tmp = a[i]; a[i] = a[j]; a[j] = tmp;
        i++; j--;
      }
    }
    if(k <= j)
      hi = j;
    else if(k >= i)
      lo = i;
    else
      return a[k];
  }
  return a[k];
}

// Median as numpy takes it: the mean of the two middle values when n is even
static double median(double *a,int n){
  if(n <= 0)
    return 0;
  double const upper = select_kth(a,n,n/2);
  if(n & 1)
    return upper;
  // after selection every value below index n/2 is <= upper; the lower middle is their max
  double lower = a[0];
  for(int i=1; i < n/2; i++)
    if(a[i] > lower)
      lower = a[i];
  return (lower + upper) / 2;
}

// One template's profile, peak and SNR over the linear span
static void do_template(struct tick_state *t,int k,int64_t window_start,struct tick_second *out){
  float const * const c = t->tmpl[k];
  int const L = t->tlen[k];
  float complex * const prof = out->profile[k];
  float const *base = (float const *)(t->lin + (t->win_off - t->span_begin));
  for(int i=0; i < t->nprofile; i++){
    float const *x = base + 2 * i * t->step;
    float re = 0, im = 0;
    for(int j=0; j < L; j++){
      re += x[2*j] * c[j];
      im += x[2*j+1] * c[j];
    }
    prof[i] = CMPLXF(re,im);
  }
  int kmax = 0;
  double pmax = -1;
  for(int i=0; i < t->nprofile; i++){
    double const p = (double)crealf(prof[i]) * crealf(prof[i]) + (double)cimagf(prof[i]) * cimagf(prof[i]);
    t->power[i] = p;
    if(p > pmax){
      pmax = p;
      kmax = i;
    }
  }
  double delta = 0;
  if(kmax > 0 && kmax < t->nprofile - 1){
    double const a = sqrt(t->power[kmax-1]), b = sqrt(t->power[kmax]), d = sqrt(t->power[kmax+1]);
    double const den = a - 2*b + d;
    if(den < 0)
      delta = 0.5 * (a - d) / den;
  }
  out->peak_index[k] = (double)window_start + (kmax + delta) * t->step;
  // Median power of the lags more than SNR_EXCLUDE_S from the peak
  int n = 0;
  for(int i=0; i < t->nprofile; i++){
    if(abs(i - kmax) * t->step > t->excl)
      t->power[n++] = t->power[i];
  }
  double const noise = median(t->power,n);
  if(noise > 0)
    out->peak_snr_db[k] = 10 * log10(pmax / noise);
  else
    out->peak_snr_db[k] = pmax > 0 ? 300.0 : 0.0; // no noise floor at all (synthetic input)
}

static void do_markers(struct tick_state const *t,struct tick_second *out){
  float complex const *x = t->lin - t->span_begin; // nominal second
  int const N = t->marker_n;
  for(int m=0; m < TICK_MARKERS; m++){
    double const w = 2*M_PI*Marker_hz[m] / t->cfg.samprate;
    double const coeff = 2 * cos(w);
    double complex s1 = 0, s2 = 0;
    for(int n=0; n < N; n++){
      double complex const s0 = x[n] + coeff * s1 - s2;
      s2 = s1;
      s1 = s0;
    }
    // |X(+w)| = |s[N-1] - e^{-iw} s[N-2]|, |X(-w)| = |s[N-1] - e^{+iw} s[N-2]|
    double complex const e = cexp(-I*w);
    double const plus = cabs(s1 - e * s2), minus = cabs(s1 - conj(e) * s2);
    out->marker_energy[m] = (plus * plus + minus * minus) / ((double)N * N);
  }
}

static void do_envelope(struct tick_state const *t,struct tick_second *out){
  float const * const m = t->mag - t->span_begin; // nominal second
  _Static_assert(TICK_ENVELOPE == ENV_RATE,"envelope is one second at 1 kHz");
  double env[TICK_ENVELOPE];
  double big = 0;
  for(int i=0; i < TICK_ENVELOPE; i++){
    float const *x = m + i * t->decim - t->env_half;
    double acc = 0;
    for(int k=0; k < t->env_ntaps; k++)
      acc += t->env_taps[k] * x[k];
    env[i] = acc;
    if(fabs(acc) > big)
      big = fabs(acc);
  }
  if(big == 0){
    out->envelope_scale = 0;
    memset(out->envelope,0,sizeof(out->envelope));
    return;
  }
  float const scale = (float)(big / 32767);
  out->envelope_scale = scale;
  for(int i=0; i < TICK_ENVELOPE; i++){
    long q = lrint(env[i] / scale);
    if(q > 32767)
      q = 32767;
    else if(q < -32767)
      q = -32767;
    out->envelope[i] = (int16_t)q;
  }
}

static void make_report(struct tick_state *t,int64_t second,struct tick_second *out){
  int64_t const begin = second + t->span_begin;
  out->second_index = second;
  out->window_start_index = second + t->win_off;
  out->step = t->step;
  out->nprofile = t->nprofile;
  out->complete = begin >= t->valid_from && begin >= t->end - t->ringsize
    && second + t->span_end <= t->end;
  if(!out->complete){
    memset(out->profile,0,sizeof(out->profile));
    for(int k=0; k < TICK_NTEMPLATES; k++){
      out->peak_index[k] = NAN;
      out->peak_snr_db[k] = 0;
    }
    memset(out->marker_energy,0,sizeof(out->marker_energy));
    memset(out->envelope,0,sizeof(out->envelope));
    out->envelope_scale = 0;
    return;
  }
  // Linear copy of the span, and its magnitude
  int p = ring_pos(t,begin);
  for(int i=0; i < t->span_len; i++){
    float complex const s = t->ring[p];
    t->lin[i] = s;
    t->mag[i] = sqrtf(crealf(s) * crealf(s) + cimagf(s) * cimagf(s));
    if(++p == t->ringsize)
      p = 0;
  }
  for(int k=0; k < TICK_NTEMPLATES; k++){
    if(t->cfg.templates & (1u << k)){
      do_template(t,k,out->window_start_index,out);
    } else {
      memset(out->profile[k],0,sizeof(out->profile[k]));
      out->peak_index[k] = NAN;
      out->peak_snr_db[k] = 0;
    }
  }
  do_markers(t,out);
  do_envelope(t,out);
}

int tick_process(struct tick_state *t,float complex const *s,int n,int64_t first_index,
                 int64_t second_index,struct tick_second *out){
  if(t == NULL || (s == NULL && n > 0) || n < 0 || out == NULL)
    return -1;
  int ret = 0;
  if(!t->have_end){
    t->have_end = true;
    t->valid_from = first_index;
  } else if(first_index != t->end){
    if(first_index < t->end){
      // Counter went backwards: everything queued refers to the old count
      t->npending = 0;
      t->have_announced = false;
    }
    t->valid_from = first_index;
    ret = -1;
  }
  t->end = first_index;
  // Into the ring
  int p = ring_pos(t,first_index);
  for(int i=0; i < n; i++){
    t->ring[p] = s[i];
    if(++p == t->ringsize)
      p = 0;
  }
  t->end = first_index + n;
  // Queue a newly named second (at least half a second after the last one)
  if(!t->have_announced || second_index >= t->last_announced + t->cfg.samprate / 2){
    if(t->npending == PENDING){
      memmove(t->pending,t->pending + 1,(PENDING - 1) * sizeof(t->pending[0]));
      t->npending--;
    }
    t->pending[t->npending++] = second_index;
    t->have_announced = true;
    t->last_announced = second_index;
  }
  if(ret < 0)
    return ret;
  // The oldest queued second, once the stream has passed its span
  if(t->npending > 0 && t->pending[0] + t->span_end <= t->end){
    make_report(t,t->pending[0],out);
    t->npending--;
    memmove(t->pending,t->pending + 1,t->npending * sizeof(t->pending[0]));
    return 1;
  }
  return 0;
}

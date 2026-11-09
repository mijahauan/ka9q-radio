// Time-station tick detector: per-second complex correlation profiles of the
// WWV, WWVH and BPM tick templates, tone-marker energies and the AM envelope
// Pure: no radiod headers, no allocation after tick_create(), no threads of its own
// Copyright 2026, HamSCI
#ifndef _TICK_H
#define _TICK_H 1
#include <stdint.h>
#include <stdbool.h>
#include <complex.h>

enum tick_template { TICK_WWV, TICK_WWVH, TICK_BPM, TICK_NTEMPLATES };

// Most profile values any accepted sample rate needs: 70 ms of lags at the
// step nearest 0.1 ms is at most 980 values (14 kHz, step 1)
#define TICK_MAX_PROFILE 1050
#define TICK_ENVELOPE 1000   // envelope values per second (1 kHz)
#define TICK_MARKERS 3       // marker tones 1000, 1200, 1500 Hz

struct tick_config {
  int samprate;        // channel output rate; a multiple of 1000, at least 4000
  unsigned templates;  // bitmask of 1<<tick_template
};

struct tick_second {                 // one report
  int64_t window_start_index;        // channel-output counter of the window's first sample
  int step;                          // samples per profile step
  int nprofile;                      // values per profile
  float complex profile[TICK_NTEMPLATES][TICK_MAX_PROFILE];
  double peak_index[TICK_NTEMPLATES];   // channel-output counter, fractional (parabolic on |profile|)
  double peak_snr_db[TICK_NTEMPLATES];
  double marker_energy[TICK_MARKERS];   // 1000, 1200, 1500 Hz over 0-800 ms
  int16_t envelope[TICK_ENVELOPE]; float envelope_scale;
  bool complete;                     // false if the window was not fully covered
};

struct tick_state;
struct tick_state *tick_create(struct tick_config const *cfg);
void tick_destroy(struct tick_state *t);
void tick_reset(struct tick_state *t);
// Feed n consecutive samples starting at channel-output counter first_index; second_index is the
// channel-output counter of the nominal second this block belongs to (the caller computes it).
// Returns 1 and fills *out when a second's report is complete, 0 otherwise, -1 on a gap.
int tick_process(struct tick_state *t, float complex const *s, int n, int64_t first_index,
                 int64_t second_index, struct tick_second *out);

char const *tick_template_name(enum tick_template k); // "wwv", "wwvh", "bpm"
// Samples past the nominal second a report needs before it can complete
int tick_span_end(struct tick_state const *t);

#endif

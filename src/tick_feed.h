// Pure glue between radiod and the tick core (tick.c): which UTC second each
// block names, where its window goes on the channel-output counter, input
// epochs, and the slot through which a pilot channel lends its time to the
// tick channels of the same radiod
// Pure: no radiod headers, no allocation, no threads of its own
// Copyright 2026, HamSCI
#ifndef _TICK_FEED_H
#define _TICK_FEED_H 1
#include <stdint.h>
#include <stdbool.h>
#include <stdatomic.h>
#include <complex.h>
#include "tick.h"
#include "pilot.h"

// Channel-output counter of UTC second utc_second_ns, given a reference pair:
// channel-output index ref_index holds UTC ref_utc_ns. Rounded to the nearest sample
int64_t tick_second_index(int64_t utc_second_ns, int64_t ref_index, int64_t ref_utc_ns, double samprate);

// The preset value "wwv", "wwv+wwvh", "wwv+wwvh+bpm" ... (any '+' list of
// wwv, wwvh, bpm, any case) as a template mask; 0 for empty, NULL or an unknown name
unsigned tick_parse_templates(char const *s);

// A reference pair for placing windows: the input-counter position whose
// content is UTC instant utc_ns as it enters radiod, before any channel filter
struct tick_ref {
  double input_pos;
  int64_t utc_ns;
  uint64_t epoch;          // the input epoch input_pos belongs to
};

// The latest locked pilot report of this radiod, as a reference pair. One
// writer at a time (writers spin on wlock for a few stores), any number of
// readers, none of which ever waits: a seqlock whose fields are all atomic
// The reference a pilot report gives: its raw edge, which runs late by the
// pilot channel's filter delay, less that delay. False (nothing written) for
// a report that is not valid and locked
bool tick_ref_from_pilot(struct pilot_report const *r, double in_samprate, struct tick_ref *ref);

struct tick_pilot_slot {
  atomic_flag wlock;
  _Atomic uint64_t seq;    // odd while being written; 0 = never written
  _Atomic int64_t input_pos;
  _Atomic int64_t utc_ns;
  _Atomic uint64_t epoch;
};
#define TICK_PILOT_SLOT_INIT { .wlock = ATOMIC_FLAG_INIT }
void tick_slot_publish(struct tick_pilot_slot *s, struct tick_ref const *r);
// False when nothing was ever published, or when writes kept tearing the read
bool tick_slot_read(struct tick_pilot_slot *s, struct tick_ref *r);

#define TICK_PILOT_MAX_AGE_S 120
// A pilot reference serves a block at input index now_input only from the
// same epoch and at most TICK_PILOT_MAX_AGE_S old
bool tick_ref_fresh(struct tick_ref const *r, uint64_t epoch, int64_t now_input, double in_samprate);

// Per-channel state. The channel-output counter is (input index - base) / decim,
// with base the input index of the channel's first block; it never rebases on
// an epoch change, so it stays monotonic with radiod's input counter
#define TICK_FEED_NAMED 8
struct tick_feed {
  int decim;               // input samples per output sample
  double in_samprate;
  int samprate;            // channel output rate
  double delay_out;        // channel filter group delay, output samples (window placement only)
  bool based;
  int64_t base;
  bool have_epoch;
  uint64_t epoch;          // epoch of the last block fed
  // Seconds named to the core, mirrored with tick.c's half-second rule so
  // each report finds the UTC second and the reference that named it
  bool have_named;
  int64_t last_named;
  struct { int64_t second_index, utc_ns; uint64_t epoch; bool used_pilot; } named[TICK_FEED_NAMED];
  int next_named;
};
struct tick_meta {         // what the glue knows about one report beside struct tick_second
  int64_t second_utc_ns;   // the UTC second measured, ns (GPS_TIME's scale)
  uint64_t epoch;          // epoch of the block that completed the report
  bool used_pilot;         // the second was named from a pilot reference, not the host clock
};
void tick_feed_init(struct tick_feed *f, int decim, double in_samprate, int samprate, double delay_out);
// One block of n channel samples whose first sample is input index
// block_input_index, in input epoch `epoch`; straddle: the block holds both
// sides of a sample loss. ref places the window; ref_is_pilot says where it
// came from. Returns 1 and fills *out and *meta when a report is complete,
// 0 otherwise, -1 for a gap or a block kept out of the core
int tick_feed_block(struct tick_state *t, struct tick_feed *f, float complex const *s, int n,
                    int64_t block_input_index, uint64_t epoch, bool straddle,
                    struct tick_ref const *ref, bool ref_is_pilot,
                    struct tick_second *out, struct tick_meta *meta);

#endif

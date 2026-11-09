// Pure glue between radiod and the tick core: naming each block's second,
// placing the window, input epochs, and the pilot slot. See tick_feed.h.
//
// Naming the second. Each block names to the core the first UTC second whose
// window (-10 ms) has not yet begun at the block's first sample, so a second
// is named about a second ahead and never after its window has started. The
// UTC instant of the block's first sample comes from a reference pair (struct
// tick_ref): a UTC instant and the input-counter position where its content
// enters radiod. The core keeps a second only if it falls at least half a
// second after the last one it kept (tick.c), so once named, a second keeps
// its alignment even if the reference moves by less than half a second
// before its window (a switch from the host clock to the pilot, say). The
// feed mirrors that rule, so each report carries the UTC second and the
// reference (pilot or host) that actually placed its window.
//
// Placing the window, and the filter delays. A channel's output sample k
// carries the input around base + k*decim minus the channel filter's group
// delay: an event that enters at input position X appears at output position
// (X - base)/decim + delay_out, with delay_out = (M-1)/2 output samples
// (pilot_filter_delay_s() in pilot_chan.c). The pilot's raw edge runs late by
// the pilot channel's own delay, so whoever builds a pilot reference first
// takes that delay off the edge (tick_pilot_publish in tick_chan.c) to get the
// input position where the UTC second entered; the host-clock reference pairs
// gps_time_ns() with the input sample count, which no filter has touched. Here
// the tick channel's own delay goes back on. With both delays the same, as on
// channels built alike, they cancel. This only moves the window, which needs
// millisecond accuracy (the ticks sit 0-50 ms inside a 70 ms window). Every
// position the core reports (second_index, window_start_index, peak_index)
// stays a raw channel-output counter value, exact on radiod's input counter
// through pilot_out_to_input and never adjusted by a model.
//
// Epochs. The output counter follows radiod's input counter, which counts
// every sample that arrives and never goes back, so it does not rebase on a
// new epoch; instead tick_break makes every report whose span began before the
// epoch's first block come out incomplete. A block that straddles a loss stays
// out of the core, and the next block then arrives as a gap.
//
// Copyright 2026, HamSCI
#include <ctype.h>
#include <strings.h>
#include <math.h>
#include <stdatomic.h>
#include <stddef.h>
#include <string.h>
#include "tick_feed.h"

#define NS 1000000000LL
#define WINDOW_LEAD_NS 10000000LL  // tick.c WINDOW_BEGIN_S; the core's span begins no earlier
#define SLOT_READ_TRIES 4

int64_t tick_second_index(int64_t utc_second_ns,int64_t ref_index,int64_t ref_utc_ns,double samprate){
  // The difference holds exactly in a double up to 2^53 ns (104 days)
  return ref_index + llround((double)(utc_second_ns - ref_utc_ns) * samprate / (double)NS);
}

unsigned tick_parse_templates(char const *s){
  if(s == NULL)
    return 0;
  unsigned mask = 0;
  while(true){
    size_t len = strcspn(s,"+");
    enum tick_template k;
    for(k = 0; k < TICK_NTEMPLATES; k++){
      char const *name = tick_template_name(k);
      if(strlen(name) == len && strncasecmp(s,name,len) == 0)
        break;
    }
    if(len == 0 || k == TICK_NTEMPLATES)
      return 0; // empty or unknown name
    mask |= 1u << k;
    if(s[len] == '\0')
      return mask;
    s += len + 1;
  }
}

void tick_slot_publish(struct tick_pilot_slot *s,struct tick_ref const *r){
  while(atomic_flag_test_and_set_explicit(&s->wlock,memory_order_acquire))
    ;
  uint64_t const q = atomic_load_explicit(&s->seq,memory_order_relaxed);
  atomic_store_explicit(&s->seq,q + 1,memory_order_relaxed); // odd: writing
  atomic_thread_fence(memory_order_release);
  atomic_store_explicit(&s->input_pos,llround(r->input_pos),memory_order_relaxed);
  atomic_store_explicit(&s->utc_ns,r->utc_ns,memory_order_relaxed);
  atomic_store_explicit(&s->epoch,r->epoch,memory_order_relaxed);
  atomic_store_explicit(&s->seq,q + 2,memory_order_release);
  atomic_flag_clear_explicit(&s->wlock,memory_order_release);
}

bool tick_slot_read(struct tick_pilot_slot *s,struct tick_ref *r){
  for(int i=0; i < SLOT_READ_TRIES; i++){
    uint64_t const q1 = atomic_load_explicit(&s->seq,memory_order_acquire);
    if(q1 == 0)
      return false; // never written
    if(q1 & 1)
      continue;     // a writer is mid-way
    int64_t const pos = atomic_load_explicit(&s->input_pos,memory_order_relaxed);
    int64_t const utc = atomic_load_explicit(&s->utc_ns,memory_order_relaxed);
    uint64_t const epoch = atomic_load_explicit(&s->epoch,memory_order_relaxed);
    atomic_thread_fence(memory_order_acquire);
    if(atomic_load_explicit(&s->seq,memory_order_relaxed) == q1){
      r->input_pos = (double)pos;
      r->utc_ns = utc;
      r->epoch = epoch;
      return true;
    }
  }
  return false; // torn every time: the caller falls back on the host clock this block
}

bool tick_ref_fresh(struct tick_ref const *r,uint64_t epoch,int64_t now_input,double in_samprate){
  if(r == NULL || r->epoch != epoch || !(in_samprate > 0))
    return false;
  return ((double)now_input - r->input_pos) <= TICK_PILOT_MAX_AGE_S * in_samprate;
}

void tick_feed_init(struct tick_feed *f,int decim,double in_samprate,int samprate,double delay_out){
  memset(f,0,sizeof(*f));
  f->decim = decim;
  f->in_samprate = in_samprate;
  f->samprate = samprate;
  f->delay_out = delay_out;
}

static int64_t floor_div(int64_t a,int64_t b){
  int64_t q = a / b;
  if((a % b != 0) && ((a < 0) != (b < 0)))
    q--;
  return q;
}

int tick_feed_block(struct tick_state *t,struct tick_feed *f,float complex const *s,int n,
                    int64_t block_input_index,uint64_t epoch,bool straddle,
                    struct tick_ref const *ref,bool ref_is_pilot,
                    struct tick_second *out,struct tick_meta *meta){
  if(t == NULL || f == NULL || ref == NULL || out == NULL || meta == NULL || f->decim < 1 || f->samprate <= 0)
    return -1;
  if(!f->based){
    f->based = true;
    f->base = block_input_index;
  }
  int64_t d = block_input_index - f->base;
  if(d < 0 || d % f->decim != 0){
    // Off the output grid (radiod's integer block ratio rules it out): start over here
    tick_reset(t);
    f->base = block_input_index;
    f->have_named = false;
    f->next_named = 0;
    memset(f->named,0,sizeof(f->named));
    d = 0;
  }
  int64_t const first = d / f->decim; // channel-output counter of this block's first sample
  if(f->have_epoch && epoch != f->epoch)
    tick_break(t,first);
  f->have_epoch = true;
  f->epoch = epoch;
  if(straddle)
    return -1;

  // Where the reference's UTC instant appears on this channel's output counter
  double const ref_out = (ref->input_pos - (double)f->base) / f->decim + f->delay_out;
  double const per_sample_ns = (double)NS / f->samprate;
  int64_t const t_first = ref->utc_ns + llround(((double)first - ref_out) * per_sample_ns);
  int64_t const second_ns = (floor_div(t_first + WINDOW_LEAD_NS - 1,NS) + 1) * NS; // ceil
  // An integer reference index for tick_second_index, its UTC carried with it
  int64_t const ref_index = (int64_t)floor(ref_out);
  int64_t const ref_utc = ref->utc_ns - llround((ref_out - (double)ref_index) * per_sample_ns);
  int64_t const second_index = tick_second_index(second_ns,ref_index,ref_utc,f->samprate);

  // Mirror tick.c's rule for keeping a named second
  if(!f->have_named || second_index >= f->last_named + f->samprate / 2){
    f->named[f->next_named].second_index = second_index;
    f->named[f->next_named].utc_ns = second_ns;
    f->named[f->next_named].used_pilot = ref_is_pilot;
    f->next_named = (f->next_named + 1) % TICK_FEED_NAMED;
    f->have_named = true;
    f->last_named = second_index;
  }
  int const r = tick_process(t,s,n,first,second_index,out);
  if(r != 1)
    return r;
  meta->second_utc_ns = 0;
  meta->used_pilot = false;
  meta->epoch = epoch;
  for(int i=0; i < TICK_FEED_NAMED; i++){
    if(f->named[i].second_index == out->second_index && f->named[i].utc_ns != 0){
      meta->second_utc_ns = f->named[i].utc_ns;
      meta->used_pilot = f->named[i].used_pilot;
      break;
    }
  }
  return 1;
}

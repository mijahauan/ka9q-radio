// radiod glue for the time-station tick detector: the "ticks" channel option,
// the pilot slot, and feeding each block to the core through tick_feed.c.
//
// Threads. tick_on_block runs in the channel's demod thread, once per fresh
// block, after pilot_on_block has recorded the anchor. The core does its work
// there: about 0.3 ms of CPU once a second per channel for all three templates
// at 24 kHz (tick_test's test_report_cost), far inside a 20 ms block, so it
// needs no worker. chan->ticks is written only by the demod thread, which also
// sends the channel's status (response() in radio.c), so status sees it consistent.
//
// The pilot slot. radiod runs one front end per process (Frontend in radio.c),
// so one slot per process is one per front end. A pilot channel's demod
// thread publishes each locked report into it; tick channels' demod threads
// read it lock-free (a seqlock, tick_feed.c), never waiting on the writer. A
// reference from another epoch or older than 120 s does not serve: the block
// then names its second from the host clock (gps_time_ns() paired with the
// input sample count), and chan->ticks.used_pilot records which one served.
//
// Copyright 2026, HamSCI
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "misc.h"
#include "radio.h"
#include "timing.h"
#include "tick.h"
#include "tick_feed.h"
#include "pilot_chan.h"
#include "tick_chan.h"

static struct tick_pilot_slot Pilot_slot = TICK_PILOT_SLOT_INIT;

void tick_pilot_publish(struct pilot_report const * const r,double const in_samprate){
  struct tick_ref ref;
  if(tick_ref_from_pilot(r,in_samprate,&ref)) // valid and locked only
    tick_slot_publish(&Pilot_slot,&ref);
}

// Set up the detector on the first block. False = can't run; the caller turns
// ticks off for this channel (until the next demod restart)
static bool tick_chan_start(chan_t * const chan){
  struct frontend const * const frontend = chan->frontend;
  if(chan->filter2.blocking != 0){
    fprintf(stderr,"%s ticks: off; needs filter2 disabled (baseband must come from the channel filter)\n",chan->name);
    return false;
  }
  int const samprate = chan->output.samprate;
  double const ratio = samprate > 0 ? frontend->samprate / samprate : 0;
  if(!(ratio >= 1) || ratio != floor(ratio) || ratio > 1e9){
    fprintf(stderr,"%s ticks: off; input rate %.1f Hz / output rate %d Hz is not an integer\n",
	    chan->name,frontend->samprate,samprate);
    return false;
  }
  struct tick_config const cfg = { .samprate = samprate, .templates = chan->ticks.templates };
  struct tick_state * const det = tick_create(&cfg);
  if(det == NULL){
    fprintf(stderr,"%s ticks: off; rate %d Hz must be a multiple of 1000 Hz, at least 4000\n",chan->name,samprate);
    return false;
  }
  struct tick_feed * const feed = calloc(1,sizeof(*feed));
  struct tick_second * const last = calloc(1,sizeof(*last)); // ~27 KB: off chan_t
  if(feed == NULL || last == NULL){
    free(feed);
    free(last);
    tick_destroy(det);
    return false;
  }
  double const delay_out = pilot_filter_delay_s(chan) * samprate;
  tick_feed_init(feed,(int)ratio,frontend->samprate,samprate,delay_out);
  chan->ticks.det = det;
  chan->ticks.feed = feed;
  chan->ticks.last = last;
  chan->ticks.decim = (int)ratio;
  chan->ticks.last_valid = false;
  fprintf(stderr,"%s ticks: on; templates%s%s%s, decimation %d, filter delay %.1f samples\n",chan->name,
	  (cfg.templates & (1u << TICK_WWV)) ? " wwv" : "",
	  (cfg.templates & (1u << TICK_WWVH)) ? " wwvh" : "",
	  (cfg.templates & (1u << TICK_BPM)) ? " bpm" : "",
	  (int)ratio,delay_out);
  return true;
}

void tick_chan_stop(chan_t * const chan){
  tick_destroy(chan->ticks.det);
  free(chan->ticks.feed);
  free(chan->ticks.last);
  chan->ticks.det = NULL;
  chan->ticks.feed = NULL;
  chan->ticks.last = NULL;
  chan->ticks.last_valid = false;
}

void tick_on_block(chan_t * const chan,bool const straddle){
  if(chan->ticks.templates == 0 || chan->baseband == NULL || chan->sampcount <= 0)
    return;
  if(chan->ticks.det == NULL && !tick_chan_start(chan)){
    chan->ticks.templates = 0;
    return;
  }
  struct frontend const * const frontend = chan->frontend;
  int64_t const index = chan->anchor.input_sample;
  uint64_t const epoch = chan->anchor.epoch;

  struct tick_ref ref;
  bool const pilot = tick_slot_read(&Pilot_slot,&ref) && tick_ref_fresh(&ref,epoch,index,frontend->samprate);
  if(!pilot){
    // The host clock now, paired with the input samples counted at that moment
    // (as pilot_on_block names the pilot's second). USB and callback latency
    // make it name UTC slightly late, so windows open early by a few ms,
    // inside the window's 10 ms lead
    ref.input_pos = (double)frontend->samples;
    ref.utc_ns = gps_time_ns();
    ref.epoch = epoch;
  }
  struct tick_meta meta;
  int const r = tick_feed_block(chan->ticks.det,chan->ticks.feed,chan->baseband,chan->sampcount,
				index,epoch,straddle,&ref,pilot,chan->ticks.last,&meta);
  if(r != 1)
    return;
  chan->ticks.last_valid = true;
  chan->ticks.last_seq++;
  chan->ticks.last_epoch = meta.epoch;
  chan->ticks.last_second_utc_ns = meta.second_utc_ns;
  chan->ticks.used_pilot = meta.used_pilot;
  chan->ticks.base_input = chan->ticks.feed->base;
}

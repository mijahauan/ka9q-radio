// radiod glue for the time-station tick detector (tick.c, tick_feed.c)
// Copyright 2026, HamSCI
#ifndef _TICK_CHAN_H
#define _TICK_CHAN_H 1
#include <stdbool.h>
#include "radio.h"

// Called from pilot_on_block on every fresh block, after chan->anchor holds
// the block's input index and epoch; straddle: the block spans a sample loss.
// With ticks.templates set, feeds the detector and collects reports into chan->ticks
void tick_on_block(chan_t *chan, bool straddle);
// Called from pilot_chan_stop as the demod exits: frees the detector
void tick_chan_stop(chan_t *chan);
// Called by a pilot channel with each locked report: lends its time to the
// tick channels of this radiod
void tick_pilot_publish(struct pilot_report const *r, double in_samprate);

#endif

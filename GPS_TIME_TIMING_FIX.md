# GPS_TIME Timing Fix for Precision Metrology

## Problem

If I understand your code (and I may well not), it creates the `GPS_TIME` field in radiod status packets at status packet generation time, not at the time corresponding to `INPUT_SAMPLES`. This introduces a systematic timing offset equal to the FFT pipeline latency (~70ms at typical settings) that will also vary between RTP streams.

## Root Cause

In `radio_status.c`, `GPS_TIME` was set by calling `gps_time_ns()` at the moment the status
packet was being encoded. However, `INPUT_SAMPLES` (from `chan->filter.out.sample_index`) 
refers to samples that arrived earlier and have passed through the FFT pipeline.

The mismatch: GPS_TIME = "now", but INPUT_SAMPLES = "samples that arrived earlier"

## Solution

Capture GPS time at the moment samples arrive (in the USB callback), store it with the
corresponding sample count, and compute GPS_TIME for any sample index using:

```
GPS_TIME = gps_time_snapshot + (sample_index - samples_at_snapshot) * BILLION / sample_rate
```

The snapshot is updated every second to track any clock drift.

### Changes

1. **radio.h**: Added `gps_time_snapshot` and `samples_at_snapshot` fields to `struct frontend`
2. **All frontend drivers**: Update snapshot every second in sample callback (at sample arrival time)
   - rx888.c, airspy.c, airspyhf.c, bladerf.c, fobos.c, funcube.c, hackrf.c, hydrasdr.c, rtlsdr.c, sdrplay.c, sig_gen.c
3. **radio_status.c**: Report `GPS_TIME` and `RTP_TIMESNAP` from snapshot (identical for all channels)
4. **linear.c**: Initialize RTP timestamp from `filter.out.sample_index / decimation` (accounts for pipeline delay)

## Why Periodic Updates?

A single capture at startup would accumulate drift over time. By updating the snapshot
every second:
- Maximum drift is bounded within ~1 second worth of samples
- Clock drift (NTP adjustments, etc.) is tracked
- No need to restart radiod to refresh timing

## Impact

- Eliminates ~70ms systematic offset in GPS_TIME/INPUT_SAMPLES relationship
- Enables sub-millisecond timing accuracy for metrology applications
- Enables coherent timing between multiple radiod instances sharing the same GPSDO and GPS+PPS reference
- No change to existing radiod behavior for non-timing-critical applications

## Testing

After applying this patch and restarting radiod:
1. Time signal arrivals (WWV, WWVH, CHU) should match expected propagation delays
2. WWV should arrive before WWVH (Fort Collins is closer than Hawaii)
3. Timing errors should be <50ms (ionospheric variation), not ~500ms (pipeline latency)

## Future HF PPS injection

Shouldn't affect anything in that regard.

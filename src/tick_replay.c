// tick-replay FILE.cf32 RATE MASK SECOND0_INDEX
// Feed a little-endian interleaved float32 I/Q capture through the tick
// detector in 20 ms batches (first_index = file sample index), streaming,
// never holding the file. SECOND0_INDEX is the file sample index of the first
// whole second; each batch names the second holding its first sample
// (SECOND0_INDEX before that). MASK: bitmask of templates (1 wwv, 2 wwvh, 4 bpm).
// Prints one JSON line per report (make_golden.py's key names)
// Copyright 2026, HamSCI
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sysexits.h>
#include "tick.h"

static void num(char const *key,double v){
  if(isfinite(v))
    printf("\"%s\": %.17g, ",key,v);
  else
    printf("\"%s\": null, ",key); // NaN is not JSON
}

int main(int argc,char *argv[]){
  if(argc < 5){
    fprintf(stderr,"Usage: %s FILE.cf32 RATE MASK SECOND0_INDEX\n",argv[0]);
    exit(EX_USAGE);
  }
  struct tick_config const cfg = {
    .samprate = atoi(argv[2]),
    .templates = (unsigned)strtoul(argv[3],NULL,0),
  };
  int64_t const second0 = strtoll(argv[4],NULL,0);
  FILE *fp = fopen(argv[1],"rb");
  if(fp == NULL){
    perror(argv[1]);
    exit(1);
  }
  struct tick_state *t = tick_create(&cfg);
  if(t == NULL){
    fprintf(stderr,"bad config: rate %d mask %u\n",cfg.samprate,cfg.templates);
    exit(1);
  }
  int const batch = cfg.samprate / 50;
  float complex *buf = malloc(batch * sizeof(*buf)); // float complex = two little-endian float32 on our targets
  static struct tick_second r;
  if(buf == NULL)
    exit(1);
  int64_t index = 0;
  size_t n;
  while((n = fread(buf,sizeof(buf[0]),batch,fp)) > 0){
    int64_t second = second0;
    if(index > second0)
      second = second0 + ((index - second0) / cfg.samprate) * cfg.samprate;
    int const ret = tick_process(t,buf,(int)n,index,second,&r);
    if(ret < 0){
      fprintf(stderr,"tick_process error at index %lld\n",(long long)index);
      exit(1);
    }
    if(ret == 1){
      printf("{\"second_index\": %lld, \"window_start_index\": %lld, \"step\": %d, \"nprofile\": %d, ",
             (long long)r.second_index,(long long)r.window_start_index,r.step,r.nprofile);
      for(int k=0; k < TICK_NTEMPLATES; k++){
        if(!(cfg.templates & (1u << k)))
          continue;
        char const *name = tick_template_name(k);
        printf("\"profile_%s\": [",name);
        for(int i=0; i < r.nprofile; i++)
          printf("%s[%.9g, %.9g]",i ? ", " : "",crealf(r.profile[k][i]),cimagf(r.profile[k][i]));
        printf("], ");
        char key[64];
        snprintf(key,sizeof(key),"peak_index_%s",name);
        num(key,r.peak_index[k]);
        snprintf(key,sizeof(key),"peak_snr_db_%s",name);
        num(key,r.peak_snr_db[k]);
      }
      printf("\"marker_energy\": [%.17g, %.17g, %.17g], ",
             r.marker_energy[0],r.marker_energy[1],r.marker_energy[2]);
      printf("\"envelope\": [");
      for(int i=0; i < TICK_ENVELOPE_LEN; i++)
        printf("%s%d",i ? ", " : "",r.envelope[i]);
      printf("], \"envelope_scale\": %.9g, \"complete\": %s}\n",
             r.envelope_scale,r.complete ? "true" : "false");
      fflush(stdout);
    }
    index += n;
  }
  fclose(fp);
  free(buf);
  tick_destroy(t);
  exit(0);
}

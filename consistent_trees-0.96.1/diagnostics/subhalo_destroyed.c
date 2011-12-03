#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include "../src/check_syscalls.h"
#include "stringparse.h"
#include <math.h>

#define EXTRA_HALO_INFO float macc, vacc, mpeak, vpeak;
#include "../read_tree/read_tree.h"
#include "../src/litehash.h"

#define BOX_SIZE 250

struct halo_list all_halos = {0};
struct halo_tree halo_tree = {0};

void _calc_subhalo_stats(struct halo *h, FILE *output)
{
  struct halo *parent = h->uparent;
  float pos[3], r=0;
  int64_t k;

  if (!parent->mvir || !parent->rvir || !parent->vmax || !parent->rs) return;
  if (h->mmp) return;

  for (k=0; k<3; k++) {
    pos[k] = h->pos[k] - parent->pos[k];
    if (pos[k] > (BOX_SIZE/2)) pos[k] = BOX_SIZE - pos[k];
    if (pos[k] < -(BOX_SIZE/2)) pos[k] = BOX_SIZE + pos[k];
  }
  for (k=0; k<3; k++) r += pos[k]*pos[k];
  r = sqrt(r);

  fprintf(output, "%"PRId64" %"PRId64" %.3e %.3e %.3g %.3g\n",
	  h->id, parent->id, h->mvir/parent->mvir,
	  parent->mvir,  r*1e3/parent->rvir, h->mpeak / parent->mpeak);
}

void calc_subhalo_stats()
{
  int64_t i;
  fprintf(stdout, "#ID UPID Mass/Pmass Pmass R/Prvir Mpeak/PMpeak\n");
  for (i=0; i<all_halos.num_halos; i++) {
    if (!all_halos.halos[i].uparent) continue;
    _calc_subhalo_stats(all_halos.halos+i, stdout);
  }
}


void read_hlist(char *filename) {
  int64_t n;
  FILE *input;
  struct halo h = {0};
  int64_t desc_pid, descid, desc_scale;
  char buffer[1024];
  struct litehash *lh;

  SHORT_PARSETYPE;
  #define NUM_INPUTS 27
  enum short_parsetype stypes[NUM_INPUTS] = 
    { F, D64, F, D64, D64,    //  #scale id desc_scale desc_id num_prog
      D64, D64, D64, D64,       //   pid upid desc_pid phantom 
      F, F, F, F, F,    //mvir orig_mvir rvir rs vrms 
      D64, F, F,          //mmp? scale_of_last_MM vmax 
      F, F, F, F, F, F, //x y z vx vy vz
      F, F, F, F};      //macc mpeak vacc vpeak
  enum parsetype types[NUM_INPUTS];
  void *data[NUM_INPUTS] = {&(h.scale), &(h.id), &(desc_scale),
                            &(descid), &(h.num_prog), &(h.pid),
			    &(h.upid), &(desc_pid), &(h.phantom), 
			    &(h.mvir), &(h.orig_mvir), &(h.rvir), &(h.rs), &(h.vrms),
			    &(h.mmp), &(h.scale_of_last_MM), &(h.vmax),
			    &(h.pos[0]), &(h.pos[1]), &(h.pos[2]), 
			    &(h.vel[0]), &(h.vel[1]), &(h.vel[2]),
			    &(h.macc), &(h.mpeak), &(h.vacc), &(h.vpeak)};
  

  for (n=0; n<NUM_INPUTS; n++) types[n] = stypes[n];
  input = check_fopen(filename, "r");
  while (fgets(buffer, 1024, input)) {
    if (buffer[0] == '#') continue;
    n = stringparse(buffer, data, (enum parsetype *)types, NUM_INPUTS);
    if (n<NUM_INPUTS) continue;
    if (!(all_halos.num_halos%3000))
      all_halos.halos = check_realloc(all_halos.halos, sizeof(struct halo)*(all_halos.num_halos+3000), "Allocating Halos.");
   
    all_halos.halos[all_halos.num_halos] = h;
    all_halos.num_halos++;
  }
  fclose(input);
  
  all_halos.halos = check_realloc(all_halos.halos, sizeof(struct halo)*all_halos.num_halos, "Allocating Halos.");

  lh = new_litehash(4);
  for (n=0; n<all_halos.num_halos; n++)
    lh_setval(lh, &(all_halos.halos[n].id), all_halos.halos + n);

  for (n=0; n<all_halos.num_halos; n++) {
    all_halos.halos[n].uparent = lh_getval(lh, &(all_halos.halos[n].upid));
    all_halos.halos[n].parent = lh_getval(lh, &(all_halos.halos[n].pid));
  }
  free_litehash(lh);
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    printf("Usage: %s hlist\n", argv[0]);
    exit(1);
  }
  read_hlist(argv[1]);
  calc_subhalo_stats();
  return 0;
}



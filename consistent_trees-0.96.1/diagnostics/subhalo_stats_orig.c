#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include "../src/check_syscalls.h"
#include "stringparse.h"
#include <math.h>

#include "../read_tree/read_tree.h"
#include "../src/litehash.h"

#define BOX_SIZE 250

struct halo_list all_halos = {0};
struct halo_tree halo_tree = {0};

double dot(float *x1, float *x2) {
  int64_t i;
  double s=0;
  for (i=0; i<3; i++) s+=x1[i]*x2[i];
  return s;
}

void decompose_v(float *pos, float *vel, float *vr, float *vtheta) {
  float vnorm;
  float r = sqrt(dot(pos, pos));
  *vr = *vtheta = 0;
  if (!r) return;
  *vr = dot(vel, pos) / r;

  vnorm = dot(vel, vel);
  *vtheta = vnorm*vnorm - (*vr)*(*vr);
  if (*vtheta > 0) *vtheta = sqrt(*vtheta);
  else *vtheta = 0;
}

void _calc_subhalo_stats(struct halo *h, FILE *output)
{
  struct halo *parent = h->parent;
  float pos[3], vel[3], r=0,  vr, vtheta;
  int64_t k;

  if (!parent->mvir || !parent->rvir || !parent->vmax || !parent->rs) return;

  for (k=0; k<3; k++) {
    pos[k] = h->pos[k] - parent->pos[k];
    if (pos[k] > (BOX_SIZE/2)) pos[k] = BOX_SIZE - pos[k];
    if (pos[k] < -(BOX_SIZE/2)) pos[k] = BOX_SIZE + pos[k];
  }
  for (k=0; k<3; k++) vel[k] = h->vel[k] - parent->vel[k];
  for (k=0; k<3; k++) r += pos[k]*pos[k];
  r = sqrt(r);
  //theta = calc_angle(pos, parent->J);
  decompose_v(pos, vel, &vr, &vtheta);
  fprintf(output, "%"PRId64" %"PRId64" %.3e %.3e %.3f %.3f %.3f %.3f %.1f "
	  "%.3f %.3f %.3f\n", h->id, parent->id, h->mvir/parent->mvir,
	  parent->mvir, h->vmax/parent->vmax, parent->vmax, h->rvir/parent->rvir, 
	  parent->rvir/1e3, parent->rvir / parent->rs, r*1e3/parent->rvir,
	  vr/parent->vmax, vtheta/parent->vmax);
}

void calc_subhalo_stats()
{
  int64_t i;
  fprintf(stdout, "#ID UPID Mass/Pmass Pmass Vmax/Pvmax Pvmax Rvir/Prvir Prvir Pc R/Prvir Vr/Pvmax Vtheta/Pvmax\n");
  for (i=0; i<all_halos.num_halos; i++) {
    if (!all_halos.halos[i].parent) continue;
    _calc_subhalo_stats(all_halos.halos+i, stdout);
  }
}


#define GROUP_LIST all_halos.halos
#define RADIUS rvir
#define FAST3TREE_TYPE struct halo
#include "fast3tree.c"
#define parent pid
#include "parents.c"
#undef parent

void read_hlist(char *filename) {
  int64_t n;
  FILE *input;
  struct halo h = {0};
  int64_t descid, np;
  char buffer[1024];
  struct litehash *lh;

  SHORT_PARSETYPE;
  #define NUM_INPUTS 14
  enum short_parsetype stypes[NUM_INPUTS] = 
    { D64, D64, F, F, F,    //  #id desc_id mvir vmax vrms
      F, F, D64, F,       //  Rvir Rs Np x y z vx vy vz 
      F, F, F, F, F,    //mvir orig_mvir rvir rs vrms 
    };
  enum parsetype types[NUM_INPUTS];
  void *data[NUM_INPUTS] = {&(h.id),
                            &(descid),			    
			    &(h.mvir), &(h.vmax), &(h.vrms), &(h.rvir), &(h.rs), 
			    &(np), 
			    &(h.pos[0]), &(h.pos[1]), &(h.pos[2]), 
			    &(h.vel[0]), &(h.vel[1]), &(h.vel[2])};
  

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

  find_parents(all_halos.num_halos);

  lh = new_litehash(4);
  for (n=0; n<all_halos.num_halos; n++)
    lh_setval(lh, &(all_halos.halos[n].id), all_halos.halos + n);

  for (n=0; n<all_halos.num_halos; n++) {
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



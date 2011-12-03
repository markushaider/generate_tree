#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include "gravitational_consistency.h"
#include "masses.h"
#include "check_syscalls.h"
#include "stringparse.h"

extern float min_mvir, max_mvir;
extern float box_size;

void gzip_file(char *filename) {
  char buffer[1024];
  snprintf(buffer, 1024, "gzip -f \"%s\"", filename);
  system(buffer);
}

void print_halo(FILE *o, struct tree_halo *th) {
  int64_t mmp;
  int64_t flags;
  if (!th) {
    fprintf(o, "#ID DescID Mvir Vmax Vrms Rvir Rs Np X Y Z VX VY VZ Jx Jy Jz spin Phantom MMP Suspicious? PID UPID Tracked Tracked_Single_MMP Num_MMP_Phantoms Original_ID\n");
     fprintf(o, "#Units: Masses in Msun / h\n"
            "#Units: Positions in Mpc / h (comoving)\n"
            "#Units: Velocities in km / s (physical)\n"
            "#Units: Angular Momenta in (Msun/h) * (Mpc/h) * km/s (physical)\n"
            "#Units: Radii in kpc / h (comoving)\n");
    return;
  }

  mmp = (th->flags & MMP_FLAG) ? 1 : 0;
  flags = 0;
  if (th->flags & SUSPICIOUS_LINK_FLAG) flags = 1;
  if (th->flags & UNPHYSICAL_LINK_FLAG) flags = 2;
  if (th->flags & MERGER_FLUCTUATION_FLAG) flags = 3;
  flags += (th->flags & (TOO_MANY_PHANT_FLAG | SHORT_TRACK_FLAG | MISTRACKED_SUB_FLAG));
  fprintf(o, "%"PRId64" %"PRId64" %.3e %.2f %.2f %.3f %.3f %"PRId64" %.5f %.5f %.5f %.2f %.2f %.2f %.3e %.3e %.3e %.5f %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64"\n",
	  th->id, th->descid, th->mvir, th->vmax, th->vrms, th->rvir, th->rs, th->np,
	  th->pos[0], th->pos[1], th->pos[2], th->vel[0], th->vel[1],
	  th->vel[2], th->J[0], th->J[1], th->J[2], th->spin, th->phantom, mmp, flags, th->pid, th->upid, th->tracked,
	  th->tracked_single_mmp, th->num_mmp_phantoms, th->orig_id);
}

void common_print_halos(char *filename, struct tree_halo *h, int64_t num_halos, int64_t gzip_halos)
{
  int64_t i;
  FILE *o;

  o = check_fopen(filename, "w");
  print_halo(o, NULL);
  for (i=0; i<num_halos; i++) print_halo(o, &(h[i]));
  fclose(o);
  if (gzip_halos) gzip_file(filename);
}

int really_beyond_mmp_ratio(struct tree_halo h1, struct tree_halo h2) {
  float m1, m2, v1, v2;
  if (h1.mvir > h2.mvir) { m1 = h1.mvir; m2 = h2.mvir; }
  else { m1 = h2.mvir; m2 = h1.mvir; }
  if (m1*0.3*MIN_MMP_MASS_RATIO > m2) return 1;

  if (h1.vmax > h2.vmax) { v1 = h1.vmax; v2 = h2.vmax; }
  else { v2 = h1.vmax; v1 = h2.vmax; }
  if (v1*0.5*MIN_MMP_VMAX_RATIO > v2) return 1;

  return 0;
}

int beyond_mmp_ratio(struct tree_halo h1, struct tree_halo h2) {
  float m1, m2, v1, v2;
  if (h1.mvir > h2.mvir) { m1 = h1.mvir; m2 = h2.mvir; }
  else { m1 = h2.mvir; m2 = h1.mvir; }
  if (m1*MIN_MMP_MASS_RATIO > m2) return 1;

  if (h1.vmax > h2.vmax) { v1 = h1.vmax; v2 = h2.vmax; }
  else { v2 = h1.vmax; v1 = h2.vmax; }
  if (v1*MIN_MMP_VMAX_RATIO > v2) return 1;

  return 0;
}

void build_id_conv_list(struct halo_stash *h)
{
  int64_t n, ids;
  if (h->num_halos < 1) return;

  h->min_id = h->max_id = h->halos[0].id;
  for (n=0; n<h->num_halos; n++) {
    if (h->max_id < h->halos[n].id) h->max_id = h->halos[n].id;
    if (h->min_id > h->halos[n].id) h->min_id = h->halos[n].id;
  }

  ids = (1 + h->max_id - h->min_id);
  h->id_conv = (int64_t *)check_realloc(h->id_conv, sizeof(int64_t)*ids, "halo index conversion");
  for (n=0; n<ids; n++) h->id_conv[n]=-1;

  for (n=0; n<h->num_halos; n++)
    h->id_conv[h->halos[n].id - h->min_id] = n;

  assert(h->min_id >= 0); //Otherwise, logic problems occur in rest of code.
}


void add_halo(struct halo_stash *h, struct tree_halo d)
{
  if (((h->num_halos)%1000)==0) {
    int64_t array_size = h->num_halos+1000;
    h->halos = (struct tree_halo *)
      check_realloc(h->halos,sizeof(struct tree_halo)*array_size, "Adding new halo.");
  }

  h->halos[h->num_halos] = d;
  h->num_halos++;
}

void load_halos(char *filename, struct halo_stash *h, float scale, int dead)
{
  FILE *input;
  char buffer[1024];
  struct tree_halo d = {0};
  int64_t n;
  int64_t mmp, flags, phantom;
  float delta_mvir = delta_vir(scale);
  float mean_density = 2.77519737e11*Om; //(Msun/h) / (comoving Mpc/h)^3
  float vir_density = delta_mvir*mean_density;
  SHORT_PARSETYPE;
#define NUM_INPUTS 27
  enum short_parsetype stypes[NUM_INPUTS] = 
    { D64, D64, F, F, F, F, F, D64, F, F, F, F, F, F, F, F, F, F, D64, D64, D64, D64, D64, D64, D64, D64, D64 };
  enum parsetype types[NUM_INPUTS];
  void *data[NUM_INPUTS] = {&(d.id),
			    &(d.descid), &(d.mvir), &(d.vmax), &(d.vrms),
			    &(d.rvir), &(d.rs), &(d.np), &(d.pos[0]),
			    &(d.pos[1]), &(d.pos[2]), &(d.vel[0]), 
			    &(d.vel[1]), &(d.vel[2]), &(d.J[0]), &(d.J[1]), &(d.J[2]),
			    &(d.spin),
			    &(phantom), &(mmp), &(flags), &(d.pid), &(d.upid),
			    &(d.tracked), &(d.tracked_single_mmp), 
			    &(d.num_mmp_phantoms), &(d.orig_id)};

  gen_ff_cache();
  d.num_prog = d.flags = d.phantom = 0;
  d.mmp_id = -1;
  for (n=0; n<NUM_INPUTS; n++) types[n] = stypes[n];
  
  input = check_fopen(filename, "r");
  while (fgets(buffer, 1024, input)) {
    if (buffer[0] == '#') continue;
    mmp = flags = 0;
    n = stringparse(buffer, data, (enum parsetype *)types, NUM_INPUTS);
    
    if (n < 14) continue;
    if (n < 18) { d.J[0] = d.J[1] = d.J[2] = d.spin = 0; }
    if (n < 19) d.phantom = 0;
    else d.phantom = phantom;
    if (n < 22) d.pid = -1;
    if (n < 23) d.upid = -1;
    if (n < 24) d.tracked = 0;
    if (n < 25) d.tracked_single_mmp = 0;
    if (n < 26) d.num_mmp_phantoms = 0;
    if (n < 27) d.orig_id = -1;

    d.mvir = fabs(d.mvir);
    if (!(d.mvir>0)) continue;
    if (!(d.rvir > 0))
      d.rvir = cbrt(d.mvir / (4.0*M_PI*vir_density/3.0)) * 1000.0;
    if (!(d.rs > 0)) d.rs = d.rvir / concentration(d.mvir, scale);

    d.flags = 0;
    if (mmp) d.flags |= MMP_FLAG;
    if ((flags & 3)==1) d.flags |= SUSPICIOUS_LINK_FLAG;
    if ((flags & 3)==2) d.flags |= UNPHYSICAL_LINK_FLAG;
    if (dead) d.flags |= DEAD_HALO_FLAG;
    add_halo(h, d);

    if (min_mvir==0) min_mvir = d.mvir;
    if (d.mvir < min_mvir) min_mvir = d.mvir;
    if (d.mvir > max_mvir) max_mvir = d.mvir;
    for (n=0; n<3; n++) if (d.pos[n] > box_size) box_size = d.pos[n];
  }
  fclose(input);

  box_size = (int) (box_size + 0.5);
}


void read_outputs(float **output_scales, int64_t **outputs, int64_t *num_outputs) {
  FILE *input;
  int64_t n;
  char buffer[1024];
  
  input = check_fopen(SCALEFILE, "r");
  while (fgets(buffer, 1024, input)) {
    if (((*num_outputs)%1000 == 0)) {
      *output_scales = (float *)check_realloc(*output_scales, sizeof(float)*((*num_outputs)+1000), "output scales");
      *outputs = (int64_t *)check_realloc(*outputs, sizeof(int64_t)*((*num_outputs)+1000), "output numbers");
    }
    n = sscanf(buffer, "%"SCNd64" %f", &((*outputs)[*num_outputs]), &((*output_scales)[*num_outputs]));
    if (n<2) continue;
    *num_outputs = (*num_outputs) + 1;
  }
  fclose(input);

  *output_scales = (float *)check_realloc(*output_scales, sizeof(float)*((*num_outputs)), "output scales");
  *outputs = (int64_t *)check_realloc(*outputs, sizeof(int64_t)*((*num_outputs)), "output numbers");
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>

#define EXTRA_HALO_INFO struct halo *vpeak_halo, *mpeak_halo, *acc_halo; int64_t sub_ok;
#include "../read_tree/read_tree.h"
#include "../read_tree/read_tree.c"
#include "../read_tree/check_syscalls.h"
#include "../src/universe_time.h"

float box_size = 250;
void process_tree(void);
float calc_dist(struct halo *h1, struct halo *h2);

int main(int argc, char **argv) {
  int64_t i;
  if (argc<3) {
    printf("Usage: %s box_size tree.dat ...\n", argv[0]);
    exit(1);
  }

  box_size = atof(argv[1]);
  init_time_table(0.27, 0.7);

  printf("#Scale Id Mnow Vnow Mpeak Vpeak Macc Vacc Rmerged(kpc/h) RHost(kpc/h) MHost VHost Rmerged/RHost Mistracked?\n");

  for (i=2; i<argc; i++) {
    read_tree(argv[i]);
    process_tree();
    delete_tree();
  }
  return 0;
}

void calc_mass_vmax_peak(struct halo *h) {
  if (h->mpeak_halo) return;
  if (h->prog) {
    calc_mass_vmax_peak(h->prog);
    
    if (h->vmax > h->prog->vpeak_halo->vmax) h->vpeak_halo = h;
    else h->vpeak_halo = h->prog->vpeak_halo;

    if (h->mvir > h->prog->mpeak_halo->mvir) h->mpeak_halo = h;
    else h->mpeak_halo = h->prog->mpeak_halo;

    if (h->upid>-1) { // If we are a subhalo
      h->acc_halo = h->prog->acc_halo;
    } else {
      h->acc_halo = h;
    }
  } else {
    h->acc_halo = h->mpeak_halo = h->vpeak_halo = h;    
  }

  if (h->acc_halo->pid < 0) h->sub_ok = 1;
  else h->sub_ok = 0;
}

void process_tree(void) {
  float rmerged, ratio;
  int64_t i;
  struct halo *h, *host;
  for (i=0; i<all_halos.num_halos; i++) {
    h = &(all_halos.halos[i]);
    if (h->pid < 0) continue;
    if (h->mmp) continue;
    if (!h->parent) continue;
    if (!h->desc) continue;
    calc_mass_vmax_peak(h);

    host = h->desc->prog;
    rmerged = calc_dist(h, host)*1e3;
    if (h->rvir>0) { ratio = rmerged / host->rvir; }
    else ratio = 0;
    
    //printf("#Scale Id Mnow Vnow Mpeak Vpeak Macc Vacc Rmerged(kpc/h) RHost(kpc/h) MHost VHost Rmerged/RHost Mistracked?\n");
    printf("%.5f %"PRId64" %.3e %.2f %.3e %.2f %.3e %.2f %f %f %.3e %.2f %f %"PRId64"\n", 
	   h->scale, h->id, h->mvir, h->vmax, h->mpeak_halo->mvir, 
	   h->vpeak_halo->vmax, h->acc_halo->mvir, h->acc_halo->vmax,
	   rmerged, host->rvir, host->mvir, host->vmax, ratio, (1-h->sub_ok));
  }
}

float calc_dist(struct halo *h1, struct halo *h2) {
  float sum = 0, dt, av_a, cmvng_mpc_per_kms, dx;
  int64_t i;
  struct halo h3;
  memcpy(&h3, h2, sizeof(struct halo));
  h2 = &h3;
  if (h2->scale != h1->scale) {
    dt = scale_to_years(h2->scale) - scale_to_years(h1->scale);
    av_a = (h1->scale+h2->scale)/2.0;
    cmvng_mpc_per_kms = 1.022e-6/av_a*dt/1.0e6;
    for (i=0; i<3; i++)
      h2->pos[i] -= cmvng_mpc_per_kms*h2->vel[i];
  }

  for (i=0; i<3; i++) {
    dx = fabs(h1->pos[i] - h2->pos[i]);
    if (dx > box_size / 2.0) dx = box_size - dx;
    sum += dx * dx;
  }
  return sqrt(sum);
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>

#define EXTRA_HALO_INFO struct halo *last_mm_halo, *acc_halo, *last_majm_halo, *last_millim_halo; int64_t sub_ok;
#include "../read_tree/read_tree.h"
#include "../read_tree/read_tree.c"
#include "../read_tree/check_syscalls.h"
#include "../src/universe_time.h"

float box_size = 250;
void process_tree(void);
float calc_radius(struct halo *h, struct halo *p);
void calc_host_mvir_rvir(struct halo *h, struct halo *p, float *mvir_acc, float *rvir_acc);

int main(int argc, char **argv) {
  int64_t i;
  if (argc<3) {
    printf("Usage: %s box_size tree.dat ...\n", argv[0]);
    exit(1);
  }

  box_size = atof(argv[1]);
  init_time_table(0.27, 0.7);

  printf("#Scale Id Mvir_acc(Host) Rvir_acc(Host) Mnow Vnow Mlastmm Vlastmm Rlastmm Rlastmm/Rvir Vacc Mistracked? RlastMajM Rlastmillim a_acc a_lastmm a_lastmajm a_lastmillim\n");

  for (i=2; i<argc; i++) {
    read_tree(argv[i]);
    process_tree();
    delete_tree();
  }
  return 0;
}

void calc_last_minor_merger(struct halo *h) {
  struct halo *prog;
  if (h->prog) {
    calc_last_minor_merger(h->prog);
    
    prog = h->prog->next_coprog;
    h->last_mm_halo = h->prog->last_mm_halo;
    h->last_majm_halo = h->prog->last_majm_halo;
    h->last_millim_halo = h->prog->last_millim_halo;
    while (prog && !(prog == h->prog)) {
      calc_last_minor_merger(prog);
      if (prog->acc_halo->mvir > 0.01 * h->prog->mvir) {
	h->last_millim_halo = h;
	if (prog->acc_halo->mvir > 0.1 * h->prog->mvir) {
	  h->last_mm_halo = h;
	  if (prog->acc_halo->mvir > 0.3 * h->prog->mvir) {
	    h->last_majm_halo = h;
	  }
	}
      }
      prog = prog->next_coprog;
    }

    if (h->upid>-1) { // If we are a subhalo
      h->acc_halo = h->prog->acc_halo;
    } else {
      h->acc_halo = h;
    }
  } else {
    h->last_majm_halo = h->last_millim_halo = h->acc_halo = h->last_mm_halo = h;    
  }

  if (h->acc_halo->pid < 0) h->sub_ok = 1;
  else h->sub_ok = 0;
}

void process_tree(void) {
  float rlast_mm, //rmpeak, 
    mvir_acc, rvir_acc, ratio, rlast_majm, rlast_millim;
  int64_t i;
  struct halo *h;
  for (i=0; i<all_halos.num_halos; i++) {
    h = &(all_halos.halos[i]);
    if (h->pid < 0) continue;
    calc_last_minor_merger(h);
    if (!h->parent) continue;
    
    rlast_mm = calc_radius(h->last_mm_halo, h->parent);
    rlast_majm = calc_radius(h->last_majm_halo, h->parent);
    rlast_millim = calc_radius(h->last_millim_halo, h->parent);
    //rmpeak = calc_radius(h->mpeak_halo, h->parent);
    calc_host_mvir_rvir(h->acc_halo, h->parent, &mvir_acc, &rvir_acc);

    ratio = (rvir_acc > 0) ? (rlast_mm / rvir_acc) : 0.0;
    printf("%.5f %"PRId64" %.3e %.4f %.3e %.2f %.3e %.2f %f %f %f %"PRId64" %.3f %.3f %.4f %.4f %.4f %.4f\n", h->scale, h->id, mvir_acc, rvir_acc, h->mvir, h->vmax, h->last_mm_halo->mvir, h->last_mm_halo->vmax, rlast_mm, ratio, h->acc_halo->vmax, (1-h->sub_ok), rlast_majm, rlast_millim, h->acc_halo->scale, h->last_mm_halo->scale, h->last_majm_halo->scale, h->last_millim_halo->scale);
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

struct halo* find_earliest_parent(struct halo *h) {
  struct halo *p = h;
  while (p->prog) p = p->prog;
  return p;
}

struct halo *find_parent_at_scale(float scale, struct halo *h) {
  struct halo *p = h;
  while (p->scale > scale) {
    p = p->prog;
    if (!p) return 0;
  }
  return p;
}

float calc_radius(struct halo *h, struct halo *p) {
  struct halo *p2 = find_parent_at_scale(h->scale, p);
  if (!p2) p2 = find_earliest_parent(p);
  if (!p2) return 0;
  return calc_dist(h, p2);
}

void calc_host_mvir_rvir(struct halo *h, struct halo *p, float *mvir_acc, float *rvir_acc) {
 struct halo *p2 = find_parent_at_scale(h->scale, p);
  if (!p2) p2 = find_earliest_parent(p);
  if (!p2) { *mvir_acc = *rvir_acc = 0; return; }
  *mvir_acc = p2->mvir;
  *rvir_acc = p2->rvir/1.0e3;
}


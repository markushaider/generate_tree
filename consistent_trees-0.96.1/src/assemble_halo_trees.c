#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include "halo_io.h"
#include "cached_io.h"
#include "litehash.h"
#include "grav_config.h"
#include "gravitational_consistency.h"
#include "gravitational_consistency_vars.h"
#include "stringparse.h"
#include "assemble_halo_trees.h"
#include "check_syscalls.h"
#include "version.h"

FILE **tree_outputs = NULL;
float *output_scales=NULL;
int64_t *outputs=NULL, num_outputs=0;
int64_t *num_trees = NULL;
struct merger_halo *extra_halos = NULL;
struct cached_io **tree_inputs = NULL;
struct merger_halo *halos = NULL;
struct litehash *lh = NULL;
int64_t num_halos=0;
float box_size, min_mvir, max_mvir;


int main(int argc, char **argv)
{
  int64_t i, j, k;
  int64_t tree_header_offset;
  char buffer[1024];
  struct cached_io *input;
  FILE *output;

  if (argc > 1) grav_config(argv[1], 0);
  else { 
    fprintf(stderr, "Consistent Trees, Version %s\n", TREE_VERSION);
    fprintf(stderr, "(C) 2011, Peter Behroozi.  See the LICENSE file for redistribution details.\n");
    fprintf(stderr, "Usage: %s options.cfg\n", argv[0]); exit(1); }
  read_outputs(&output_scales, &outputs, &num_outputs);

  tree_header_offset = create_headers();
  tree_inputs = check_realloc(NULL, sizeof(struct cached_io *)*num_outputs,
			      "Allocating tree inputs.");
  extra_halos = check_realloc(NULL, sizeof(struct merger_halo)*num_outputs,
			      "Allocating extra halo slots.");
  for (i=0; i<num_outputs; i++) extra_halos[i].id = -1;
  for (i=0; i<num_outputs; i++) {
    snprintf(buffer, 1024, "%s/really_consistent_%"PRId64".list", 
	     OUTBASE, outputs[i]);
    tree_inputs[i] = cfopen(buffer, 2000000);
  }

  input = tree_inputs[num_outputs-1];

  while (cfgets(input, buffer, 1024)) {
    struct merger_halo halo;
    if (buffer[0] == '#') continue;
    read_halo_from_line(&halo, buffer);
    halo.scale = output_scales[num_outputs-1];
    halo.desc_scale = 0;

    int64_t id = halo_pos_to_tree(&halo);
    num_trees[id]++;

    num_halos = 0;
    if (lh) { free_litehash2(lh); }
    lh = new_litehash(8);
    halos = add_to_array(halos, &num_halos, sizeof(struct merger_halo), &halo);

    build_tree(halo.id, num_outputs-2);
    print_tree_halos(id);
  }

  //Print num_tree counts
  for (i=0; i<BOX_DIVISIONS; i++) {
    for (j=0; j<BOX_DIVISIONS; j++) {
      for (k=0; k<BOX_DIVISIONS; k++) {
	int64_t id = i*BOX_DIVISIONS*BOX_DIVISIONS + j*BOX_DIVISIONS + k;
	output = tree_outputs[id];
	fseek(output, tree_header_offset, SEEK_SET);
	fprintf(output, "%-12"PRId64"\n", num_trees[id]);
	fclose(output);
      }
    }
  }

  return 0;
}

void print_tree_halos(int64_t file_id) {
  int64_t i;
  
  for (i=0; i<num_halos; i++)
    lh_setval2(lh, &(halos[i].scale), &(halos[i].id), &(halos[i]));

  calc_desc_pid_and_num_progs();
  calc_mmp();
  smooth_mass();
  conserve_mass();
  calc_mmp();
  calc_last_mm();

  fprintf(tree_outputs[file_id], "#tree %"PRId64"\n", halos[0].id);
  for (i=0; i<num_halos; i++)
    print_tree_halo(halos + i, tree_outputs[file_id]);
}


void calc_desc_pid_and_num_progs(void) {
  int64_t i;
  for (i=0; i<num_halos; i++) {
    double desc_scale = halos[i].desc_scale;
    int64_t desc_id = halos[i].descid;
    struct merger_halo *desc = lh_getval2(lh, &desc_scale, &desc_id);
    if (!desc) continue;
    halos[i].desc = desc;
    desc->num_prog++;
    halos[i].desc_pid = desc->pid;
  }
}

void calc_last_mm(void) {
  int64_t i;
  for (i=0; i<num_halos; i++) halos[i].last_mm = 0;
  for (i=num_halos-1; i>=0; i--) {
    if ((halos[i].last_mm==0) && (halos[i].mmp_halo) && (halos[i].mmp_halo->last_mm > 0)) {
      halos[i].last_mm = halos[i].mmp_halo->last_mm;
    }
    if (halos[i].desc && (halos[i].desc->mmp_halo->id != halos[i].id) &&
	halos[i].mvir > MAJOR_MERGER*halos[i].desc->mmp_halo->mvir) {
      halos[i].desc->last_mm = halos[i].desc->scale;
    }
  }
}


void calc_mmp(void) {
  int64_t i;
  for (i=0; i<num_halos; i++) {
    struct merger_halo *th = halos + i;
    if (!th->desc) continue;
    if (!th->desc->mmp_halo) th->desc->mmp_halo = th;
    else if (th->desc->mmp_halo->mvir < th->mvir)
      th->desc->mmp_halo = th;
  }

  for (i=0; i<num_halos; i++) {
    if (!halos[i].desc) { halos[i].mmp = 1; continue; }
    halos[i].mmp = (halos[i].desc->mmp_halo->id == halos[i].id) ? 1 : 0;
  }
}

void smooth_mass(void) {
  int64_t i;
  for (i=0; i<num_halos; i++) {
    struct merger_halo *th = halos+i;
    if (!th->desc) continue;
    if (th->desc && th->desc->mmp_halo && (th->desc->mmp_halo->id == th->id)) {
      th->next_mass = th->desc->mvir;
      th->desc->prev_mass = th->mvir;
    }
  }
  for (i=0; i<num_halos; i++) {
    struct merger_halo *th = halos+i;
    if (th->next_mass && th->prev_mass)
      th->mvir = (th->next_mass + th->prev_mass + th->mvir) / 3.0;
  }
}

int double_sort(const void *a, const void *b) {
  double c = *((double *)a);
  double d = *((double *)b);
  if (c<d) return -1;
  if (c>d) return 1;
  return 0;
}

void conserve_mass(void) {
  int64_t i, j, *ids;
  double scale, *scales = lh_keylist(lh);
  qsort(scales, lh->elems, sizeof(double), double_sort);
  for (i=0; i<num_halos; i++) halos[i].incoming_mass = 0;

  for (i=0; i<lh->elems; i++) {
    scale = scales[i];
    struct litehash *idlh = lh_getval(lh, &scale);
    assert(idlh);
    ids = lh_keylist(idlh);
    for (j=0; j<idlh->elems; j++) {
      struct merger_halo *th = lh_getval(idlh, ids+j);
      assert(th);
      if ((th->upid < 0) && th->incoming_mass
	  && (th->incoming_mass > th->mvir)) {
	th->mvir = th->incoming_mass;
      }
      if (!th->desc) continue;
	    
      if (th->pid < 0) th->desc->incoming_mass += th->mvir;
      if ((th->pid < 0) && (th->desc->pid > -1) && (th->desc->id != th->desc->pid)) {
	struct merger_halo *parent = lh_getval2(lh, &(th->desc->scale), &(th->desc->pid));
	if (parent) parent->incoming_mass += th->mvir;
      }
    }
    free(ids);
  }
  free(scales);
}

int64_t create_headers(void) {
  int64_t i, j, k, id, offset = 0;
  char buffer[1024];
  FILE *output;
  int64_t num_output_trees = BOX_DIVISIONS*BOX_DIVISIONS*BOX_DIVISIONS;
  tree_outputs = check_realloc(tree_outputs, sizeof(FILE *)*num_output_trees,
			       "Allocating outputs.");

  num_trees = check_realloc(NULL, sizeof(int64_t)*num_output_trees,
			    "Allocating tree counts.");
  memset(num_trees, 0, sizeof(int64_t)*num_output_trees);

  for (i=0; i<BOX_DIVISIONS; i++) {
    for (j=0; j<BOX_DIVISIONS; j++) {
      for (k=0; k<BOX_DIVISIONS; k++) {
	snprintf(buffer, 1024, "%s/tree_%"PRId64"_%"PRId64"_%"PRId64".dat",
		 TREE_OUTBASE, i, j, k);
	unlink(buffer);
	id = i*BOX_DIVISIONS*BOX_DIVISIONS + j*BOX_DIVISIONS + k;
	output = tree_outputs[id] = check_fopen(buffer, "w");
	fprintf(output,
		"#scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) mvir(9) orig_mvir(10) rvir(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26)\n"
		"#Scale: Scale factor of halo.\n"
		"#ID: ID of halo (unique across entire simulation).\n"
		"#Desc_Scale: Scale of descendant halo, if applicable.\n"
		"#Descid: ID of descendant halo, if applicable.\n"
		"#Num_prog: Number of progenitors.\n"
		"#Pid: Host halo ID (-1 if distinct halo).\n"
		"#Upid: Most massive host halo ID (only different from Pid in cases of sub-subs, or sub-sub-subs, etc.).\n"
		"#Desc_pid: Pid of descendant halo (if applicable).\n"
		"#Phantom: Nonzero for halos interpolated across timesteps.\n"
		"#Mvir: Halo mass, smoothed across accretion history; always greater than sum of halo masses of contributing progenitors (Msun/h).\n"
		"#Orig_Mvir: Original halo mass from raw halo catalogs (before smoothing or mass conservation), Msun/h.\n"
		"#Rvir: Halo radius (kpc/h comoving).\n"
		"#Rs: Scale radius (kpc/h comoving).\n"
		"#Vrms: Velocity dispersion (km/s physical).\n"
		"#mmp?: whether the halo is the most massive progenitor or not.\n"
		"#scale_of_last_MM: scale factor of the last major merger (Mass ratio > %g).\n"
		"#Vmax: Maxmimum circular velocity (km/s physical).\n"
		"#X/Y/Z: Halo position (Mpc/h comoving).\n"
		"#VX/VY/VZ: Halo velocity (km/s physical).\n"
		"#JX/JY/JZ: Halo angular momenta ((Msun/h) * (Mpc/h) * km/s (physical)).\n"
		"#Spin: Halo spin parameter.\n"
		"#Not the same definition as Brandon's tree: mmp\n"
		"#Please note: Do not make assumptions about the ordering of halos in each tree.\n"
		"#Consistent Trees Version %s\n", MAJOR_MERGER, TREE_VERSION);
	
	offset = ftello(output);
	fprintf(output, "XXXXXXXXXXXX\n"); //For the number of trees
      }
    }
  }
  return offset;
}

int64_t read_halo_from_line(struct merger_halo *halo, char *buffer) {
  SHORT_PARSETYPE;
#define NUM_INPUTS 23
  enum short_parsetype stypes[NUM_INPUTS] = 
    { D64, D64, F64, F64, F64, F64, F64, D64, F, F, F, F, F, F, F, F, F, F, D64, D64, K, D64, D64 };
  enum parsetype types[NUM_INPUTS];
  int64_t n;
  void *data[NUM_INPUTS] = {&(halo->id),
                            &(halo->descid), &(halo->mvir), &(halo->vmax), 
			    &(halo->vrms), &(halo->rvir), &(halo->rs), 
			    &(halo->np), &(halo->pos[0]), &(halo->pos[1]), 
			    &(halo->pos[2]), &(halo->vel[0]), &(halo->vel[1]), 
			    &(halo->vel[2]), &(halo->J[0]), &(halo->J[1]),
			    &(halo->J[2]), &(halo->spin), &(halo->phantom), 
			    &(halo->mmp), 
			    NULL, &(halo->pid), &(halo->upid)};

  memset(halo, 0, sizeof(struct merger_halo));
  for (n=0; n<NUM_INPUTS; n++) types[n] = stypes[n];
  n = stringparse(buffer, data, (enum parsetype *)types, NUM_INPUTS);
  halo->orig_mvir = halo->mvir;
  halo->desc_pid = -1;
  return ((n==NUM_INPUTS) ? 1 : 0);
#undef NUM_INPUTS
}

int64_t halo_pos_to_tree(struct merger_halo *halo) {
  int64_t idx[3], i;
  for (i=0; i<3; i++) idx[i] = (int64_t)(halo->pos[i]*BOX_DIVISIONS/BOX_WIDTH)
			% (int64_t)(BOX_DIVISIONS);
  return (idx[0]*BOX_DIVISIONS*BOX_DIVISIONS + idx[1]*BOX_DIVISIONS + idx[2]);
}

void *add_to_array(void *array, int64_t *size, int64_t width, void *data) {
  if (!((*size)%1000))
    array = check_realloc(array, width*((*size)+1000), "Allocating array elements.");
  memcpy(array + (*size)*width, data, width);
  *size = (*size) + 1;
  return array;
}

void build_tree(int64_t id, int64_t inputnum) {
  int64_t i = inputnum;
  int64_t id_index = 0, res;
  double scale, desc_scale;
  struct merger_halo halo;
  char buffer[1024];
  int64_t *ids=NULL, num_ids = 0;
  int64_t *new_ids=NULL, num_new_ids = 0;

  ids = add_to_array(ids, &num_ids, sizeof(int64_t), &id_index);

  while (i >= 0) {
    struct cached_io *input = tree_inputs[i];
    scale = output_scales[i];
    desc_scale = (i < num_outputs-1) ? output_scales[i+1] : 0;

    while (1) {
      if (extra_halos[i].id >= 0) {
	halo = extra_halos[i];
	extra_halos[i].id = -1;
      } else {
	while ((res = cfgets(input, buffer, 1024)) && (buffer[0] == '#'));
	if (!res) break;
	memset(&halo, 0, sizeof(struct merger_halo));
	if (!read_halo_from_line(&halo, buffer)) continue;
      }

      if (halo.descid < 0) {
	fprintf(stderr, "Orphaned halo found!\n");
	exit(1);
      }

      //Check ID of halo to make sure we're still in the right tree:
      for (; id_index < num_ids; id_index++) {
	if (halo.descid == halos[ids[id_index]].id) break;
      }
      if (id_index == num_ids) {
	extra_halos[i] = halo;
	break;
      }
	    
      halo.scale = scale;
      halo.desc_scale = desc_scale;
      new_ids = add_to_array(new_ids, &num_new_ids, sizeof(int64_t), &num_halos);
      halos = add_to_array(halos, &num_halos, sizeof(struct merger_halo), &halo);
    }

    free(ids);
    ids = new_ids;
    num_ids = num_new_ids;
    new_ids = NULL;
    num_new_ids = 0;
    id_index = 0;
    i--;
  }
  free(new_ids);
}


void print_tree_halo(struct merger_halo *h, FILE *output) {
  fprintf(output, " %.4f %8"PRId64" %.4f %8"PRId64" %6"PRId64" %8"PRId64" %8"PRId64" %8"PRId64" %2"PRId64" %.5e %.5e %6f %6f %6f %2"PRId64" %.4f %6f %.5f %.5f %.5f %.3f %.3f %.3f %.3e %.3e %.3e %.5f\n",
	  h->scale, h->id, h->desc_scale, h->descid, h->num_prog,
	  h->pid, h->upid, h->desc_pid, h->phantom,
	  h->mvir, h->orig_mvir, h->rvir, h->rs, h->vrms,
	  h->mmp, h->last_mm, h->vmax,
	  h->pos[0], h->pos[1], h->pos[2],
	  h->vel[0], h->vel[1], h->vel[2],
	  h->J[0], h->J[1], h->J[2], h->spin);
}


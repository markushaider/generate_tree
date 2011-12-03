#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include "gravitational_consistency.h"
#include "gravitational_consistency_vars.h"
#include "halo_io.h"
#include "grav_config.h"
#include "check_syscalls.h"
#include "masses.h"
#include "halo_io.h"
#include "resort_outputs.h"
#ifndef NO_FORK
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>
#endif /* NO_FORK */
#include "version.h"

struct halo_stash now={0}, evolved = {0}, dead = {0};
float box_size=0;
float max_mvir=0;
float min_mvir=0;
int64_t *sort_order = NULL;

int main(int argc, char **argv) {
  int64_t i, j, num_outputs=0, restart_snap = -1;
  float *output_scales=NULL;
  int64_t *outputs=NULL;
  char buffer[1024];
  int stat_loc, last_output;

  if (argc==1) {
    fprintf(stderr, "Consistent Trees, Version %s\n", TREE_VERSION);
    fprintf(stderr, "(C) 2011, Peter Behroozi.  See the LICENSE file for redistribution details.\n");
    fprintf(stderr, "Usage: %s options.cfg\n", argv[0]); exit(1);
  }
  if (argc>1) grav_config(argv[1], 1);
  read_outputs(&output_scales, &outputs, &num_outputs);
  gen_ff_cache();
  clear_halo_stash(&now);

  if (argc>2) restart_snap = i = atoi(argv[2]);
  else i = num_outputs-1;
  for (; i>=0; i--) {
    snprintf(buffer, 1024, "%s/really_consistent_%"PRId64".list", OUTBASE, outputs[i]);
    clear_halo_stash(&now);
    load_halos(buffer, &now, output_scales[i], 0);
    build_id_conv_list(&now);
    if (i==num_outputs - 1) last_output = 1;
    else if (i==restart_snap) last_output = -1;
    else last_output = 0;
    sort_order = check_realloc(sort_order, sizeof(int64_t)*now.num_halos,
			       "Allocating sort order for halos.");
    for (j=0; j<now.num_halos; j++) sort_order[j] = j;
    resort_halos(last_output);
    fork_and_print_halos(buffer, &now);

    clear_halo_stash(&dead);
    snprintf(buffer, 1024, "%s/really_dead_%"PRId64".list", OUTBASE, outputs[i]);
    load_halos(buffer, &dead, output_scales[i], 1);
    sort_order = check_realloc(sort_order, sizeof(int64_t)*dead.num_halos,
			       "Allocating sort order for halos.");
    for (j=0; j<dead.num_halos; j++) sort_order[j] = j;
    qsort(sort_order, dead.num_halos, sizeof(int64_t), sort_dead_by_mvir);
    fork_and_print_halos(buffer, &dead);
    while (wait4(-1, &stat_loc, WNOHANG, 0)>0);
  }
  while (wait4(-1, &stat_loc, 0, 0)>0);
  return 0;
}

void fork_and_print_halos(char *filename, struct halo_stash *h)
{
  char buffer[1024];
  pid_t pid;
  FILE *o;
  int64_t i;
  
  snprintf(buffer, 1024, "%s.new", filename);
  pid = fork();
  if (pid < 1) {
    o = check_fopen(buffer, "w");
    print_halo(o, NULL);
    for (i=0; i<h->num_halos; i++) print_halo(o, h->halos + sort_order[i]);      
    fclose(o);
    rename(buffer, filename);
    if (!pid) exit(0);
  }
}

inline int64_t id_to_index(struct halo_stash h, int64_t id) {
  if (id < h.min_id || id > h.max_id) return -1;
  return (h.id_conv[id-h.min_id]);
}

int sort_by_mvir(const void *a, const void *b) {
  const struct tree_halo *c = now.halos + *((int64_t *)a);
  const struct tree_halo *d = now.halos + *((int64_t *)b);
  if (c->mvir > d->mvir) return -1;
  if (d->mvir > c->mvir) return 1;
  if (c->id < d->id) return -1;
  return 1;
}

int sort_dead_by_mvir(const void *a, const void *b) {
  const struct tree_halo *c = dead.halos + *((int64_t *)a);
  const struct tree_halo *d = dead.halos + *((int64_t *)b);
  if (c->mvir > d->mvir) return -1;
  if (d->mvir > c->mvir) return 1;
  if (c->id < d->id) return -1;
  return 1;
}

int sort_by_desc(const void *a, const void *b) {
  const struct tree_halo *c = now.halos + *((int64_t *)a);
  const struct tree_halo *d = now.halos + *((int64_t *)b);
  int64_t cindex = id_to_index(evolved, c->descid);
  int64_t dindex = id_to_index(evolved, d->descid);
  if (cindex < 0 && dindex >= 0) return 1;
  if (cindex >=0 && dindex < 0) return -1;
  if (cindex < dindex) return -1;
  if (dindex < cindex) return 1;
  if (c->mvir < d->mvir) return 1;
  if (d->mvir < c->mvir) return -1;
  if (c->id < d->id) return -1;
  return 1;
}


void clear_halo_stash(struct halo_stash *h) {
  if (h->halos) h->halos = (struct tree_halo *)realloc(h->halos, 0);
  if (h->id_conv) h->id_conv = (int64_t *)realloc(h->id_conv, 0);
  h->max_id = h->min_id = h->num_halos = 0;
}

void zero_halo_stash(struct halo_stash *h) {
  h->halos = 0;
  h->id_conv = 0;
  h->max_id = h->min_id = h->num_halos = 0;
}

void resort_halos(int last_output)
{
  int64_t i;
  if (last_output > -1) {
    if (last_output)
      qsort(sort_order, now.num_halos, sizeof(int64_t), sort_by_mvir);
    else
      qsort(sort_order, now.num_halos, sizeof(int64_t), sort_by_desc);
  }
  clear_halo_stash(&evolved);
  evolved.halos = check_realloc(evolved.halos, sizeof(struct tree_halo)*now.num_halos,
				"Allocating space for sorted halos.");
  evolved.num_halos = now.num_halos;
  for (i=0; i<now.num_halos; i++) evolved.halos[i] = now.halos[sort_order[i]];
  build_id_conv_list(&evolved);
  return;
}

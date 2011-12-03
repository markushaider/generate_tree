#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include "../read_tree/read_tree.h"
#include "../read_tree/check_syscalls.h"

void add_stats(struct halo *h, int64_t sub_ok, int64_t steps_tracked, int64_t num_phantoms);
void process_tree(void);


#define MAX_BINS 200
int64_t total[MAX_BINS] = {0};
int64_t ok[MAX_BINS] = {0};

int main(int argc, char **argv) {
  int64_t i;
  FILE *output;
  if (argc<2) {
    printf("Usage: %s tree.dat ...\n", argv[0]);
    exit(1);
  }

  for (i=1; i<argc; i++) {
    read_tree(argv[i]);
    process_tree();
    delete_tree();
  }
 
  output = fopen("mistracked.dat", "w");
  for (i=0; i<MAX_BINS; i++) {
    if (!total[i]) continue;
    fprintf(output, "%f %f %"PRId64"\n", (float)i/4.0, 
	    1.0 - ((float)ok[i])/((float)total[i]), total[i]);
  }
  fclose(output);
  return 0;
}

void process_tree(void) {
  int64_t i, steps_tracked, num_phantoms, sub_ok;
  struct halo_list *hl = halo_tree.halo_lists;
  struct halo *h, *p;

  for (i=0; i<hl->num_halos; i++) {
    h = &(hl->halos[i]);
    if (h->pid < 0) continue;
    sub_ok = steps_tracked = num_phantoms = 0;
    for (p=h; p; p=p->prog) {
      steps_tracked++;
      if (p->phantom) num_phantoms++;
      if (p->pid < 0) { sub_ok = 1; break; }
    }
    add_stats(h, sub_ok, steps_tracked, num_phantoms);
  }
}

void add_stats(struct halo *h, int64_t sub_ok, int64_t steps_tracked, int64_t num_phantoms) {
  int64_t bin;
  
  if (!sub_ok)
    printf("%"PRId64" %"PRId64" %f %"PRId64" %"PRId64" %e %f %f %f %f\n", steps_tracked, num_phantoms,
	   h->scale, h->id, h->phantom, h->mvir, h->vmax,
	   h->pos[0], h->pos[1], h->pos[2]);
  if (!(h->mvir>0)) return;
  bin = 4*log10(h->mvir);
  if (bin > MAX_BINS) return;

  total[bin]++;
  if (sub_ok) ok[bin]++;
}



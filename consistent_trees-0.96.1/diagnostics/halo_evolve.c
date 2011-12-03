#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include "../minimal_halo.h"

void load_halos(char *filename, float scale);

struct min_halo *halos = NULL;
int64_t num_halos = 0;

int main(int argc, char **argv) {
  if (argc < 2) {
    printf("Usage: %s halolist1\n", argv[0]);
    exit(1);
  }

  load_halos(argv[1], 0);
  return 0;
}

void load_halos(char *filename, float scale) {
  FILE *input;
  char buffer[1024];
  struct min_halo halo = {0};
  float max_d = 0;
  int n, i;
  if (!(input = fopen(filename, "r"))) {
    printf("Couldn't open file %s!\n", filename);
    exit(2);
  }
  
  for (i=0; i<3; i++) halo.a[i] = 0;
  while (fgets(buffer, 1024, input)) {
    if (buffer[0] == '#') continue;
    n = sscanf(buffer, "%d %d %f %f %f %f %f %d %f %f %f %f %f %f %d",
	       &(halo.id),
	       &(halo.descendant), &(halo.mvir), &(halo.vmax), &(halo.vrms),
	       &(halo.rvir), &(halo.rs), &(halo.np), &(halo.pos[0]),
	       &(halo.pos[1]), &(halo.pos[2]), &(halo.vel[0]), 
	       &(halo.vel[1]), &(halo.vel[2]),
	       &(halo.phantom_id));
    if (n < 14) continue;
    if (n < 15) halo.phantom_id = 0;
    halo.mvir = fabs(halo.mvir);
    if (!(halo.mvir > 0)) continue;

    if (!(num_halos % 1000)) {
      halos = (struct min_halo *)
	realloc(halos, sizeof(struct min_halo)*(num_halos+1000));
      if (!halos) {
	printf("Out of memory trying to load halos!\n");
	exit(1);
      }
    }
    halos[num_halos] = halo;
    num_halos++;

  }
  fclose(input);
}



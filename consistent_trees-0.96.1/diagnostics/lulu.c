#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "../src/stringparse.h"

typedef struct {
  float pos[3];
  int32_t pid;
  float lum;
} halo;

#define FAST3TREE_TYPE halo
#include "../src/fast3tree.c"

halo *halos = NULL;
int64_t num_halos = 0;
struct fast3tree *tree = NULL;

void load_halos(char *filename);

#define MW_WIDTH 0.2
#define H 0.7
#define MW_LUM (-21.2 - 0.7745098) // -5 log h
#define CL_DIST (0.5*H) //Mpc
#define MW_DIST (3.1*H) //Mpc
#define SAT_DIST (0.15*H)
#define SAT_BRIGHT_CUT 2.0
#define SAT_DIM_CUT 4.0


int main(int argc, char **argv) 
{
  int64_t i, num_fg, num_fg2, num_fg3, num_sats, j, err;
  float x,y,z, r2d, r3d, max_dist;
  int box_size;
  halo *h2;
  struct fast3tree_results *res;
  if (argc < 2) {
    printf("Usage: %s halolist\n", argv[0]);
    exit(1);
  }

  res = fast3tree_results_init();
  load_halos(argv[1]);
  box_size = tree->root->max[0]+0.5;
  max_dist = ((float)box_size)/2.0;

  for (i=0; i<num_halos; i++) {
    if (halos[i].pid > -1) continue;
    if (fabs(halos[i].lum-MW_LUM) > MW_WIDTH) continue;
    fast3tree_find_sphere_periodic(tree, res, halos[i].pos, MW_DIST);
    num_fg3 = num_fg2 = err = num_fg = num_sats = 0;
    for (j=0; j<res->num_points; j++) {
      h2 = res->points[j];

#define DIST(a,b) a = h2->b - halos[i].b;	\
      if (a > max_dist) a-=box_size;		\
      else if (a < -max_dist) a+=box_size;
      
      DIST(x,pos[0]);
      DIST(y,pos[1]);
      DIST(z,pos[2]);

      r2d = x*x+y*y;
      r3d = x*x+y*y+z*z;

      if ((h2->lum < halos[i].lum) && (r2d < (CL_DIST*CL_DIST)))
	{ err=1; break; }
      if (h2->lum-SAT_BRIGHT_CUT < halos[i].lum) continue;
      if (h2->lum-SAT_DIM_CUT > halos[i].lum) continue;


      if (r2d < (SAT_DIST*SAT_DIST)) { 
	if (r3d < (SAT_DIST*SAT_DIST)) num_sats++;
	else num_fg++;
      } else if (r2d < (2*SAT_DIST*SAT_DIST)) {
	num_fg2++;
      }
      if ((r2d < (2.5*SAT_DIST*SAT_DIST)) &&
	  (r2d > (1.5*SAT_DIST*SAT_DIST)))
	num_fg3++;	  
    }
    if (!err) {
      printf("%lld %lld %lld %lld\n", num_sats, num_fg, num_fg2, num_fg3);
    }
  }
  return 0;
}

void load_halos(char *filename) {
  FILE *input;
  char buffer[1024];
  SHORT_PARSETYPE;
  halo h;
  int n;
  enum parsetype types[] = {K,K,K,K,K,
			    K,D,K,K,K,
			    K,K,K,K,K,
			    K,K,F,F,F,
			    K,K,K,K,K,
			    K,K,K,F};
  #define N NULL
  void *data[] = {N,N,N,N,N,
		  N,&(h.pid),N,N,N,
		  N,N,N,N,N,
		  N,N,&(h.pos[0]), &(h.pos[1]), &(h.pos[2]),
		  N,N,N,N,N,
		  N,N,N,&(h.lum)};
  
  input = fopen(filename, "r");
  while (fgets(buffer, 1024, input)) {
    if (buffer[0]=='#') continue;
    n = stringparse(buffer, data, types, 29);
    if (n<29) continue;

    if (!(num_halos % 1000)) {
      halos = (halo *)
        realloc(halos, sizeof(halo)*(num_halos+1000));
      if (!halos) {
        printf("Out of memory trying to load halos!\n");
        exit(1);
      }
    }

    halos[num_halos] = h;
    num_halos++;    
  }
  fclose(input);

  tree = fast3tree_init(num_halos, halos);
}

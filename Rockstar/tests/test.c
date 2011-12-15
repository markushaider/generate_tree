#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "bounds.c"

int main(void) {
  float bounds[6];
  float particle[3] = {0, 4, 0};
  int64_t in_bounds = 0;
  float newpos[3];
  float a[3], b[3];
  int64_t i;
  BOX_SIZE=1;
  PERIODIC=1;

  sscanf("-1 -1 -1 2 2 2", "%f %f %f %f %f %f", bounds, bounds+1, bounds+2,
	 bounds+3, bounds+4, bounds+5);
  for (i=0; i<1073741824*0.2; i++) {
    in_bounds+=_check_bounds(particle, newpos, bounds);
    particle[0]+=1e-9;
  }
  printf("Total: %"PRId64"\n", in_bounds);
  return 0;
}

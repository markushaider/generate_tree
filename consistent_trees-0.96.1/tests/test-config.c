#include "grav_config.h"
#include <stdio.h>
#include "gravitational_consistency_vars.h"




int main(void) {
  grav_config("default.cfg");
  printf("SIGMA_X_MIN: %f\n", SIGMA_X_MIN);
  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include "gravitational_consistency.h"
#include "read_config.h"

#define rs(x) x = config_to_string(&c,#x,x)
#define rf(x) x = config_to_real(&c,#x,x)

void grav_config(char *filename, int write_cfile) {
  struct configfile c = {0};
  char buffer[1024];
  load_config(&c, filename);

  rf(Om);
  rf(Ol);
  rf(h0);

  rs(SCALEFILE);
  rs(INBASE);
  rs(OUTBASE);
  rs(TREE_OUTBASE);

  rf(MAJOR_MERGER);
  rf(MIN_MMP_MASS_RATIO);
  rf(MIN_MMP_VMAX_RATIO);
  rf(PADDING_TIMESTEPS);

  rf(MIN_TIMESTEPS_TRACKED);
  rf(MIN_TIMESTEPS_SUB_TRACKED);
  rf(MIN_TIMESTEPS_SUB_MMP_TRACKED);
  rf(MAX_PHANTOM_FRACTION);

  rf(LAST_DITCH_SEARCH_LIMIT);
  rf(LAST_DITCH_VMAX_RATIO_1);
  rf(LAST_DITCH_VMAX_RATIO_2);

  rf(MAX_PHANTOM);
  rf(MAX_PHANTOM_SMALL);
  rf(SMALL_PARTICLE_LIMIT);
  rf(TIDAL_FORCE_LIMIT);
  rf(RECURSION_LIMIT);
  rf(METRIC_LIMIT);
  UNPHYSICAL = 3*METRIC_LIMIT + 1;
  rf(METRIC_BREAK_LIMIT);
  rf(MASS_RES_OK);

  rf(BOX_DIVISIONS);
  rf(BOX_WIDTH);
  rf(SOFTENING_LENGTH);

  if (write_cfile) {
    snprintf(buffer, 1024, "%s/grav_consistency.cfg", OUTBASE);
    write_config(c, buffer);
  }
  free_config(c);
}

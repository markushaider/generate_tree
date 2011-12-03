#ifndef _GRAVITATIONAL_CONSISTENCY_H_
#define _GRAVITATIONAL_CONSISTENCY_H_

#include "tree_halo.h"

/* STRUCTURES */
struct halo_pos {
  float pos[3];
  struct tree_halo *h;
};

struct halo_stash {
  struct tree_halo *halos;
  int64_t num_halos, min_id, max_id;
  int64_t *id_conv;
};

/* MACROS */
#ifdef NO_FORK
#define fork(x) -1
#define wait4(a,b,c,d)
typedef int64_t pid_t;
#endif /* NO_FORK */

#ifdef NO_PERIODIC
#define IF_PERIODIC if (0) 
#else
#define IF_PERIODIC if (1)
#endif /* NO_PERIODIC */

#define MMP_FLAG 1
#define DEAD_HALO_FLAG 2 //Halo to be removed from printouts
#define SUSPICIOUS_LINK_FLAG 4 //Link which is super-unlikely
#define UNPHYSICAL_LINK_FLAG 8 //Link which is technically unphysical
#define MERGER_FLUCTUATION_FLAG 16 //Halo which appears and then merges at the next timestep.
#define SHORT_TRACK_FLAG 32
#define FIND_NEW_DESC_FLAG 64
#define TOO_MANY_PHANT_FLAG 128
#define MISTRACKED_SUB_FLAG 256
#define MAJOR_MERGER_FLAG 512
#define MIGHT_BE_FAKE_FLAG 1024

/* VARIABLES */
/* See gravitational_consistency_vars.h for definitions. */
extern float Om, h0, Ol;

extern char *SCALEFILE;
extern char *INBASE;
extern char *OUTBASE;
extern char *TREE_OUTBASE;

extern float MAJOR_MERGER,
  MIN_MMP_MASS_RATIO,
  MIN_MMP_VMAX_RATIO;

extern int PADDING_TIMESTEPS, MIN_TIMESTEPS_TRACKED,
  MIN_TIMESTEPS_SUB_TRACKED, MIN_TIMESTEPS_SUB_MMP_TRACKED;

extern float MAX_PHANTOM_FRACTION;
extern float SOFTENING_LENGTH;


extern float LAST_DITCH_SEARCH_LIMIT,
  LAST_DITCH_VMAX_RATIO_1,
  LAST_DITCH_VMAX_RATIO_2;

extern int MAX_PHANTOM,
  MAX_PHANTOM_SMALL,
  SMALL_PARTICLE_LIMIT,
  RECURSION_LIMIT;

extern float BOX_WIDTH, BOX_DIVISIONS;

extern float TIDAL_FORCE_LIMIT,
  METRIC_LIMIT,
  UNPHYSICAL,
  METRIC_BREAK_LIMIT,
  MASS_RES_OK;

#endif /* _GRAVITATIONAL_CONSISTENCY_H_ */


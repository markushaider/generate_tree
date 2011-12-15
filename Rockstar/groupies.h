#ifndef GROUPIES_H
#define GROUPIES_H

#include "fof.h"
#include "halo.h"
#include <stdint.h>

#define GROWING_FLAG 1
#define DELETE_FLAG 2
#define POSSIBLE_SWAP_FLAG 4

extern struct halo *halos;
extern int64_t num_halos;
struct extra_halo_info *extra_info;

void find_subs(struct fof *f);
void calc_mass_definition(void);
void free_particle_copies(void);

//Internal functions
void add_new_halo(void);
void alloc_particle_copies(int64_t total_copies);
void norm_sd(struct fof *f, float thresh);
float find_median_r(float *rad, int64_t num_p, float frac);

void find_sd(struct particle *particles, int64_t n, double *sig_x, double *sig_v);
#endif /* GROUPIES_H */

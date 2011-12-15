#ifndef ROCKSTAR_H
#define ROCKSTAR_H

#include <stdint.h>
#include "particle.h"
#include "halo.h"
#include "fof.h"

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* From math.h */
#endif /* M_PI */

#define MIN_WORKUNIT 5000000
#define LARGE_FOF (MIN_WORKUNIT / sizeof(struct particle))

extern struct particle *p;
extern int64_t num_p;
extern int64_t num_all_fofs;

struct workunit_info {
  int64_t num_fofs, num_halos, num_particles, chunk, seed;
  float bounds[6];
};

void rockstar(float *bounds, int64_t manual_subs);
void rockstar_cleanup();
void prune_fofs(float *bounds);
void build_particle_tree(void);
void clear_particle_tree(void);
struct particle ** find_halo_sphere(struct halo *h, int64_t *num_results);
int sort_fofs(const void *a, const void *b);
int work_finished(void);
void do_workunit(struct workunit_info *w, struct fof *fofs);
void find_unfinished_workunit(struct workunit_info *w, struct fof **fofs,
			      struct particle **parts);
void integrate_finished_workunit(struct workunit_info *w, struct fof *fofs,
				 struct halo *h, struct extra_halo_info *ei,
				 struct particle *parts);

#endif /* ROCKSTAR_H */

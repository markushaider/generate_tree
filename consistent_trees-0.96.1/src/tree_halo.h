#ifndef TREE_HALO_H
#define TREE_HALO_H

#include <stdint.h>

struct tree_halo {
  int64_t id, descid, mmp_id, tidal_id, pid, upid, orig_id;
  int64_t flags, phantom;
  float tidal_force, mass_factor;
  float pos[3], vel[3], J[3], spin;
  double a[3];
  float mvir, rvir, vmax, vrms, rs;
  float pid_vmax, upid_vmax;
  int64_t np, num_prog;
  int64_t tracked, tracked_single_mmp, num_mmp_phantoms;
};

#endif /*TREE_HALO_H*/

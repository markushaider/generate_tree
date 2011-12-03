#ifndef _MERGER_HALO_H_
#define _MERGER_HALO_H_

struct merger_halo {
  int64_t id, descid, mmp, pid, upid, desc_pid;
  double scale, desc_scale;
  struct merger_halo *desc, *mmp_halo;
  int64_t phantom;
  float pos[3], vel[3], J[3], spin;
  double mvir, rvir, vmax, vrms, rs, orig_mvir;
  double next_mass, prev_mass, incoming_mass;
  float last_mm;
  float vmax_peak, vmax_acc, mvir_peak, mvir_acc;
  int64_t np, num_prog;
};

#endif /* _MERGER_HALO_H_ */

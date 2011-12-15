#ifndef _STATS_H_
#define _STATS_H_

struct fullval {
  int64_t bin;
  double val;
};

struct binstats {
  double bin_start, bin_end, bpunit;
  double *keyavg, *avg, *m2, *sd;
  double *med, *stdup, *stddn;
  struct fullval *fullvals;
  int64_t *counts;
  int64_t num_fullvals;
  int64_t fullstats, num_bins, logarithmic;
};

void clear_binstats(struct binstats *bs);
struct binstats *init_binstats(double bin_start, double bin_end, double bpunit,
			       int64_t fullstats, int64_t logarithmic);
void free_binstats(struct binstats *bs);
void add_to_binstats(struct binstats *bs, double key, double val);
void calc_median(struct binstats *bs);
void print_avgs(struct binstats *bs, FILE *output);
void print_medians(struct binstats *bs, FILE *output);

#endif /* _STATS_H_ */

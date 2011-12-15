#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include "../check_syscalls.h"
#include "stats.h"

#define ONE_SIGMA 0.341344746

void clear_binstats(struct binstats *bs) {
  if (bs->fullstats) {
    memset(bs->med, 0, sizeof(double)*(bs->num_bins+1));
    memset(bs->stdup, 0, sizeof(double)*(bs->num_bins+1));
    memset(bs->stddn, 0, sizeof(double)*(bs->num_bins+1));
    if (bs->fullvals) {
      free(bs->fullvals);
      bs->fullvals = NULL;
    }
  }
  bs->num_fullvals = 0;
  memset(bs->counts, 0, sizeof(int64_t)*(bs->num_bins+1));
  memset(bs->avg, 0, sizeof(double)*(bs->num_bins+1));
  memset(bs->sd, 0, sizeof(double)*(bs->num_bins+1));
  memset(bs->m2, 0, sizeof(double)*(bs->num_bins+1));
  memset(bs->keyavg, 0, sizeof(double)*(bs->num_bins+1));
}


struct binstats *init_binstats(double bin_start, double bin_end, double bpunit,
			       int64_t fullstats, int64_t logarithmic)
{
  struct binstats *bs = check_realloc(NULL, sizeof(struct binstats),
				      "Allocating binstats");
  bs->bin_start = bin_start;
  bs->bin_end = bin_end;
  bs->bpunit = bpunit;
  bs->num_bins = (bin_end-bin_start)*bpunit + 1.0000001;
  bs->fullstats = fullstats;
  bs->logarithmic = logarithmic;

  bs->keyavg = check_realloc(NULL, sizeof(double)*(bs->num_bins+1),
			     "Allocating keyavgs");
  bs->avg = check_realloc(NULL, sizeof(double)*(bs->num_bins+1),
			  "Allocating avgs");
  bs->m2 = check_realloc(NULL, sizeof(double)*(bs->num_bins+1),
			  "Allocating var");
  bs->sd = check_realloc(NULL, sizeof(double)*(bs->num_bins+1),
			  "Allocating var");
  if (fullstats) {
    bs->med = check_realloc(NULL, sizeof(double)*(bs->num_bins+1),
			    "Allocating medians");
    bs->stdup = check_realloc(NULL, sizeof(double)*(bs->num_bins+1),
			      "Allocating stdup");
    bs->stddn = check_realloc(NULL, sizeof(double)*(bs->num_bins+1),
			      "Allocating stddn");
  } else {
    bs->med = bs->stdup = bs->stddn = NULL;
  }

  bs->fullvals = NULL;
  bs->counts = check_realloc(NULL, sizeof(int64_t)*(bs->num_bins+1),
			  "Allocating counts");
  clear_binstats(bs);
  return bs;
}

void free_binstats(struct binstats *bs) {
  free(bs->avg);
  free(bs->keyavg);
  free(bs->sd);
  free(bs->m2);
  free(bs->counts);
  if (bs->fullvals) free(bs->fullvals);
  if (bs->med) free(bs->med);
  if (bs->stdup) free(bs->stdup);
  if (bs->stddn) free(bs->stddn);
  free(bs);
}

void add_to_binstats(struct binstats *bs, double key, double val) {
  double keyval = key, delta;
  int64_t bin;
  if (bs->logarithmic) {
    if (key <= 0) keyval = -1000;
    else keyval = log10(key);
  }
  bin = (keyval-bs->bin_start)*bs->bpunit;
  if (bin < 0 || bin > bs->num_bins) bin = bs->num_bins;
  bs->counts[bin]++;
  bs->keyavg[bin] += (keyval - bs->keyavg[bin])/(double)bs->counts[bin];
  delta = val - bs->avg[bin];
  bs->avg[bin] += delta / (double)bs->counts[bin];
  bs->m2[bin] += delta*(val - bs->avg[bin]);
  bs->sd[bin] = bs->m2[bin] /(double)bs->counts[bin];
  if (bs->fullstats) {
    if (!(bs->num_fullvals%1000))
      bs->fullvals = check_realloc(bs->fullvals, sizeof(struct fullval)*
				   (bs->num_fullvals+1000), "fullvals");
    bs->fullvals[bs->num_fullvals].bin = bin;
    bs->fullvals[bs->num_fullvals].val = val;
    bs->num_fullvals++;
  }
} 

static int sort_stats(const void *a, const void *b) {
  const struct fullval *c = a;
  const struct fullval *d = b;
  if (c->bin < d->bin) return -1;
  if (c->bin > d->bin) return 1;
  if (c->val < d->val) return -1;
  if (c->val > d->val) return 1;
  return 0;
}

void calc_median(struct binstats *bs) {
  int64_t i, v_start;
  if (!bs->fullstats) {
    fprintf(stderr, "Binstats structure needs to be set as fullstats before calculating median!\n");
    exit(1);
  }
  qsort(bs->fullvals, bs->num_fullvals, sizeof(struct fullval), sort_stats);
  v_start = 0;
  for (i=0; i<bs->num_bins; i++) {
    if (!bs->counts[i])
      bs->med[i] = bs->stdup[i] = bs->stddn[i] = 0;
    else {
      bs->med[i] = bs->fullvals[(int64_t)(v_start + bs->counts[i]*0.5)].val;
      bs->stdup[i] = bs->fullvals[(int64_t)(v_start + bs->counts[i]*(0.5+ONE_SIGMA))].val - bs->med[i];
      bs->stddn[i] = bs->med[i] - bs->fullvals[(int64_t)(v_start + bs->counts[i]*(0.5-ONE_SIGMA))].val;
    }
    v_start += bs->counts[i];
  }
}

void print_avgs(struct binstats *bs, FILE *output) {
  int64_t i;
  for (i=0; i<bs->num_bins; i++) {
    if (!bs->counts[i]) continue;
    double key = bs->keyavg[i];
    if (bs->logarithmic) key = pow(10, key);
    fprintf(output, "%g %g %g #%"PRId64"\n", key, bs->avg[i], bs->sd[i],
	    bs->counts[i]);

  }
}

void print_medians(struct binstats *bs, FILE *output) {
  int64_t i;
  calc_median(bs);
  for (i=0; i<bs->num_bins; i++) {
    if (!bs->counts[i]) continue;
    double key = bs->keyavg[i];
    if (bs->logarithmic) key = pow(10, key);
    fprintf(output, "%g %g %g %g #avg: %g; counts: %"PRId64"\n", key,
	    bs->med[i], bs->stdup[i], bs->stddn[i], bs->avg[i], bs->counts[i]);
  }
}

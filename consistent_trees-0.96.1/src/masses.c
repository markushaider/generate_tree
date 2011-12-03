#include <math.h>
#include <inttypes.h>
#include "masses.h"

#define FF_CACHE_SIZE 10001
#define HALO_PROFILE_CUTOFF 100 /* maximum value of rh/rs cached */
#define FF_CACHE_STEP (((double)HALO_PROFILE_CUTOFF)/ \
		       ((double)FF_CACHE_SIZE-1.0))
#define FF_INV_CACHE_STEP (((double)FF_CACHE_SIZE-1.0)/ \
			   ((double)HALO_PROFILE_CUTOFF))

float ff_cache[FF_CACHE_SIZE] = {0};

inline double filling_fraction(double x) {
  return (log1p(1.0/x) - 1.0/(1.0+x));
}

inline double kravtsov_f(double x) {
  return (x*x*x*filling_fraction(x));
}

double delta_vir(double a) {
  double x = 1.0/(1.0+a*a*a*Ol/Om)-1.0;
  return ((18*M_PI*M_PI + 82.0*x - 39*x*x)/(1.0+x));
}

double concentration(double mvir, double scale) {
  if (!(mvir>0)) return 0;
  double logm = log10(mvir);
  return (pow(10.0, 2.2358 - 0.10*logm)*scale);
}


/*See Appendix C of http://arxiv.org/abs/astro-ph/0203169 for mass conversions*/
double calculate_mass_factor(double mvir, double rvir, double rs)
{
  /* mvir = 4pi rho_s rvir^3 * (rs/rvir)^3 * filling_fraction(rs/rvir) */
  /* mass_factor = 4pi rho_s rs^3 = mvir / filling_fraction(rs/rvir) */
  /* Hence, m(r) = mass_factor * filling_fraction(rs/r) */
  double invc = rs / rvir;
  return (mvir / filling_fraction(invc));
}

void gen_ff_cache(void) {
  int64_t i;
  double f;
  for (i=0; i<FF_CACHE_SIZE; i++) {
    f = FF_CACHE_STEP*i;
    ff_cache[i] = i ? filling_fraction(1.0/f) : 0;
  }
}

/* Returns the filling fraction as a function of r_h / r_s */
float ff_cached(float x) {
  int64_t i;
  x *= FF_INV_CACHE_STEP;
  if (!(x>=0)) return 0;
  if (x>=FF_CACHE_SIZE-1) return ff_cache[FF_CACHE_SIZE-1];
  i = x;
  x -= i;
  return (ff_cache[i] + x*(ff_cache[i+1]-ff_cache[i]));
}


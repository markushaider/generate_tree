#include <math.h>
#include <inttypes.h>
#include "universal_constants.h"
#include "config_vars.h"

#define STEPS 1024

static double t_0 = 0;       // Time now (in Hubble time units).
static double times[STEPS+1]={0};
static double H_CONV = HUBBLE_TIME_CONVERSION/0.7;

void init_time_table(void) {
  double a = 1;
  double t = t_0;
  double dadt, dtda;
  int64_t i;
  H_CONV = HUBBLE_TIME_CONVERSION/h0;

  times[STEPS] = t_0;
  for (i=1; i<=STEPS; i++) {
    dadt = sqrt(Om * (1.0/a - (a*a)) + (a*a));
    dtda = 1.0/dadt;
    a -= 1.0/((double)STEPS);
    t -= dtda*(1.0/((double)STEPS));
    times[STEPS-i] = t;
  }
}


// Linearly interpolate between calculated values.
double scale_to_time(double scale) {
  double s = scale;
  int64_t l = (int)(s*STEPS);
  double f = s*STEPS - l;
  if (scale > 1) return ((scale-1.0)*STEPS*(times[STEPS]-times[STEPS-1]));
  if (scale < 0) return times[0];
  return (times[l]+f*(times[l+1]-times[l]));
}

double scale_to_years(double scale) {
  return (scale_to_time(scale)*H_CONV);
}

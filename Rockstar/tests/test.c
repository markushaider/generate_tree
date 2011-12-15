#include "../nfw.c"

int main(void) {
  char buffer[1024];
  double mvir, rvir, vmax, rvmax, scale, oldrs;
  while (fgets(buffer, 1024, stdin)) {
    if (sscanf(buffer, "%lf %lf %lf %lf %lf %lf", &mvir, &rvir, &vmax, &rvmax, &scale, &oldrs) < 6) continue;
    float fc = 0;
    float rs = calc_scale_radius(mvir, rvir, vmax, rvmax, scale);
    //printf("Scale radius: %f; c: %f\n", rs, rvir/rs);
    printf("%f %f %f %f %f %f %f\n", mvir, vmax, rvir/rs, rvir/oldrs, fc, rvir, vmax/cbrt(mvir));
  }
  return 0;
}

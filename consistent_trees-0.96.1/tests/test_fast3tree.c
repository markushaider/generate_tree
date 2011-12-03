#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include "distance.h"
#include "minimal_halo.h"
#include "masses.h"
#include "universe_time.h"
#include "halo_evolve.h"

#define FAST3TREE_TYPE struct min_halo 
#include "fast3tree.c"

struct min_halo *halos = NULL;
struct fast3tree *halo_tree = NULL;
int64_t num_halos = 0;
float box_size = 0; /* Automatically set; in Mpc/h */
float max_mvir = 0; /* Automatically set; in Msun */
float Om=0, Ol=0, h0=0; /* Omega_matter, Omega_lambda, h0 */

int main(int argc, char **argv) {
  float a1, a2, da, a;
  int i;

  if (argc < 6) {
    printf("Usage: %s scale1 scale2 Om h halolist1 [halolist2 ...]\n", argv[0]);
    exit(1);
  }

  /* Init cosmology and timescales */
  a1 = atof(argv[1]);
  a2 = atof(argv[2]);
  Om = atof(argv[3]);
  if (!Om) Om = 0.27;
  Ol = 1.0 - Om;
  h0 = atof(argv[4]);
  if (!h0) h0 = 0.705;
  init_cosmology(Om, Ol, h0);
  init_time_table(Om, h0);
  gen_ff_cache();

  num_halos = 0;
  max_mvir = 0;
  for (i=5; i<argc; i++) load_halos(argv[i], a1);
  halo_tree = fast3tree_init(num_halos, halos);

  da = (a2-a1)/((double)NUM_TIMESTEPS);
  do_timestep(0, a1, a1+da, 1); //First step
  return 0;
}

void load_halos(char *filename, float scale) {
  FILE *input;
  char buffer[1024];
  struct min_halo halo = {0};
  float max_d = 0;
  float delta_mvir = delta_vir(scale);
  float mean_density = 2.77519737e11*Om; //(Msun/h) / (comoving Mpc/h)^3
  float vir_density = delta_mvir*mean_density;
  int n, i;
  if (!(input = fopen(filename, "r"))) {
    printf("Couldn't open file %s!\n", filename);
    exit(2);
  }
  
  memset(halo.a, 0, sizeof(double)*3);
  while (fgets(buffer, 1024, input)) {
    if (buffer[0] == '#') continue;
    n = sscanf(buffer, "%d %d %f %f %f %f %f %d %f %f %f %f %f %f %d",
	       &(halo.id),
	       &(halo.descendant), &(halo.mvir), &(halo.vmax), &(halo.vrms),
	       &(halo.rvir), &(halo.rs), &(halo.np), &(halo.pos[0]),
	       &(halo.pos[1]), &(halo.pos[2]), &(halo.vel[0]), 
	       &(halo.vel[1]), &(halo.vel[2]),
	       &(halo.phantom_id));
    if (n < 14) continue;
    if (n < 15) halo.phantom_id = 0;
    halo.mvir = fabs(halo.mvir);
    if (!(halo.mvir > 0)) continue;
    if (!(halo.rvir > 0))
      halo.rvir = cbrt(halo.mvir / (4.0*M_PI*vir_density/3.0)) * 1000.0;
    if (!(halo.rs > 0)) halo.rs = halo.rvir / concentration(halo.mvir, scale);

    //Halo mvir is assumed to be in Msun/h
    //Note that mvir doesn't actually *have* to be the virial mass---
    //it just has to be the mass contained within whatever rvir is.
    halo.mass_factor = calculate_mass_factor(halo.mvir /h0, halo.rvir, halo.rs);
    halo.rs /= 1000.0; //Convert kpc/h to Mpc/h
    if (!(num_halos % 1000)) {
      halos = (struct min_halo *)
	realloc(halos, sizeof(struct min_halo)*(num_halos+1000));
      if (!halos) {
	printf("Out of memory trying to load halos!\n");
	exit(1);
      }
    }
    halos[num_halos] = halo;
    num_halos++;

    if (max_mvir < halo.mvir) max_mvir = halo.mvir;
    for (i=0; i<3; i++)
      if (max_d < halo.pos[i]) max_d = halo.pos[i];
  }
  fclose(input);
  box_size = (int)(max_d + 0.5); //In Mpc/h comoving coordinates.
  max_mvir /= h0; //Now in Msun
}


void set_tree_maxmin(void) {
  int i;
  for (i=0; i<3; i++) {
    halo_tree->root->min[i] = 0;
    halo_tree->root->max[i] = box_size;
  }
}

/* Assumes mvir in Msun and dt in Myr. */
inline float halo_grav_range(float mvir, float dt) {
  /* dv = dt * a = dt * (Gc*m / r^2) */
  /* So, we have dv > VEL_RESOLUTION from r=0 up to the return value: */
  return (sqrtf(fabs(mvir*Gc*dt / VEL_RESOLUTION)));
}

/* Advances simulation from a2 to a3; old a1 needed for velocity updates. */
/* special_step = 1 (first step), 0 (normal step), -1 (last step) */
void do_timestep(double a1, double a2, double a3, int special_step)
{
  int i,j,k;
  float dt = (scale_to_years(a3) - scale_to_years(a2))/1.0e6; //In Myr
  float old_dt = (scale_to_years(a2) - scale_to_years(a1))/1.0e6; //In Myr
  float range, av_a = (a2+a3)/2.0;
  float inv_rs, r,m, dx, dy, dz, r3, acc, rsoft2, inv_rs2, softening_factor;
  float dvx, dvy, dvz, vorb, coulomb, lnc, dens, rd;
  struct fast3tree_results *nearest;
  struct min_halo *h1, *h2;
  float vel_dt = dt * 1.02268944e-6 * h0 / av_a; //1 km/s to comoving Mpc/Myr/h
  float max_dist = box_size / 2.0;
  float f_h;
  int f_id;
  // 1 km/s / Myr to comoving Mpc/Myr^2/h 
  float acc_dt2 = (dt*dt / 2.0) * 1.02268944e-6 * h0 / av_a;
  float conv_const = av_a*av_a / h0 / h0;

  //Update intermediate velocities
  if (special_step != 1) { //If we're not at the first step.
    for (i=0; i<num_halos; i++) {
      h1 = &(halos[i]);
      for (j=0; j<3; j++)
	h1->vel[j] += h1->a[j]*old_dt*0.5;
      if (special_step >= 0)
	memset(h1->a, 0, sizeof(float)*3); //Clear accelerations.
    }
  }

  if (special_step<0) return; //Return if this is the final step.

  //Calculate accelerations
  if (special_step != 1) 
    fast3tree_rebuild(halo_tree, num_halos, halos);
  set_tree_maxmin();
  nearest = fast3tree_results_init();

  for (i=0; i<num_halos; i++) {
    f_h = 0;
    h1 = &(halos[i]);
    range = halo_grav_range(h1->mvir/h0, dt)*h0/av_a; //In comoving Mpc/h
    if (max_dist < range) range = max_dist*0.99;
    assert(fast3tree_find_sphere_periodic(halo_tree, nearest, h1->pos, range));

    for (j=0; j<nearest->num_points; j++) {
      r = 0;
      for (k=0; k<3; k++) {
	dx = h1->pos[k]-nearest->points[j].pos[k];
	r+=dx*dx;
      }
      r = sqrtf(r);
      if (r>f_h) { f_h = r; f_id = nearest->points[j].id; }
    }
    printf("%d %lld %f %f %f %f %f %d\n", h1->id, nearest->num_points, h1->pos[0], h1->pos[1], h1->pos[2], range, f_h, f_id);
  }
  fast3tree_results_free(nearest);
}

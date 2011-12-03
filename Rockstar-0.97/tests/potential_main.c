#define POTENTIAL_COMPARISON
#include "potential.c"

int compare_r(const void *a, const void *b) {
  const struct potential *c = a;
  const struct potential *d = b;
  if (c->r < d->r) return -1;
  if (c->r > d->r) return 1;
  return 0;
}

void calc_circ_potential(struct potential *po, int64_t num_po, float *cen)
{
  int64_t i;
  double num_bound = -0.5;
  double phi = 0;
  double potential_constant = PARTICLE_MASS;
  for (i=0; i<num_po; i++) po[i].r = sqrt(_distance2(cen, po[i].pos));
  qsort(po, num_po, sizeof(struct potential), compare_r);

  for (i=num_po-1; i>=0; i--) {
    po[i].pe2 = potential_constant*(i/po[i].r + phi);
    if (po[i].pe2 > po[i].ke)
      phi += 1.0/po[i].r;
  }
}

#define NUM_POT 2000000
int64_t num_pot = 0;
struct potential pot[NUM_POT];

int compare_x(const void *a, const void *b) {
  const struct potential *c = a;
  const struct potential *d = b;
  if (c->pos[0] < d->pos[0]) return -1;
  if (c->pos[0] > d->pos[0]) return 1;
  if (c->pos[1] < d->pos[1]) return -1;
  if (c->pos[1] > d->pos[1]) return 1;
  if (c->pos[2] < d->pos[2]) return -1;
  if (c->pos[2] > d->pos[2]) return 1;
  return 0;
}


int main(int argc, char **argv) {
  char buffer[1024];
  float f[6];
  //float cen[6] = {0};
  float cen[6] = //{26.203939, 14.303870, 6.479383, -240.270966, 103.091614, 384.215027};
    {26.203751, 14.304119, 6.480038, -225.746475, 103.026428, 409.277405};
  int64_t i;
  //PARTICLE_MASS = 1.36e8*20;
  SCALE_NOW = 1.0;
  PARTICLE_MASS = 1.36e8;
  while (fgets(buffer, 1024, stdin)) {
    sscanf(buffer, "%f %f %f %f %f %f", f, f+1, f+2, f+3, f+4, f+5);
    //for (i=0; i<6; i++) cen[i]+=f[i];
    memcpy(pot[num_pot].pos, f, sizeof(float)*6);
    pot[num_pot].flags = 0;
    num_pot++;
    if (num_pot == NUM_POT) break;
  }
  //for (i=0; i<6; i++) cen[i]/=(float)num_pot;
  //for (i=0; i<100; i++)
  compute_kinetic_energy(pot, num_pot, cen);
  compute_potential(pot, num_pot);
  calc_circ_potential(pot, num_pot, cen);
  qsort(pot, num_pot, sizeof(struct potential), compare_x);
  for (i=0; i<num_pot; i++) {
    printf("%f %f %f %g %g %g\n", pot[i].pos[0], pot[i].pos[1], pot[i].pos[2],
	   pot[i].pe, pot[i].pe2, pot[i].pe-pot[i].ke); //, pot[i].ke, 
    //(pot[i].pe > pot[i].ke) ? 1 : 0);
  }
}

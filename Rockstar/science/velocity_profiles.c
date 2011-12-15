#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>
#include "stats.h"
#include "../check_syscalls.h"
#include "../io/bgc2.h"
#include "../io/io_util.h"

GROUP_DATA_RMPVMAX *grps = NULL;
PARTICLE_DATA_PV *parts = NULL;
int64_t num_g=0, num_p=0;
int64_t num_printed = 0;

double MIN_MASS = 1e14;
double MAX_MASS = 3e14;
#define BPDEX 5
#define BIN_START (-3)
#define BIN_END 0
#define INV_BPDEX (1.0/((double)BPDEX))
#define NUM_BINS (int64_t)((BIN_END-BIN_START)*BPDEX + 1)

void load_bgc2(char *filename, struct bgc2_header *hdr,
	       GROUP_DATA_RMPVMAX **groups, int64_t *num_groups,
	       PARTICLE_DATA_PV **pdata, int64_t *num_parts);

void calc_bulkvel(float *bulkvel, PARTICLE_DATA_PV *p, int64_t length) {
  double vel[3] = {0};
  for (int64_t i=0; i<length; i++)
    for (int64_t j=0; j<3; j++) vel[j]+=p[i].vel[j];
  if (length < 1) length = 1;
  for (int64_t j=0; j<3; j++) bulkvel[j] = vel[j]/(double)length;
}

float dist(float *p1, float *p2) {
  int64_t j;
  double ds=0, dx;
  for (j=0; j<3; j++) { dx = p1[j]-p2[j]; ds+=dx*dx; }
  return sqrt(ds);
}

void bin_velocities(PARTICLE_DATA_PV *p, int64_t length, float *pos,
		    float *vel, struct binstats *bs, FILE *output) {
  int64_t i,j,bin;
  float vel_avg[3], d, dv;
  static struct binstats *this_halo[3] = {0};
  
  for (i=0; i<3; i++) {
    if (!this_halo[i])
      this_halo[i] = init_binstats(BIN_START, BIN_END, BPDEX, 0, 1);
    else clear_binstats(this_halo[i]);
  }

  for (i=0; i<length; i++) {
    d = dist(p[i].pos, pos);
    for (j=0; j<3; j++) add_to_binstats(this_halo[j], d, p[i].vel[j]);
  }
  for (bin=0; bin<NUM_BINS; bin++) {
    d = pow(10, BIN_START + (bin+0.5)*INV_BPDEX);
    if (!this_halo[0]->counts[bin]) continue;
    for (j=0; j<3; j++) vel_avg[j] = this_halo[j]->avg[bin];
    dv = dist(vel_avg, vel);
    if (output) fprintf(output, "%g %g %g %g %g %"PRId64"\n",
			d, dv, vel_avg[0]-vel[0], vel_avg[1]-vel[1],
			vel_avg[2]-vel[2], this_halo[0]->counts[bin]);
    add_to_binstats(bs, d, dv);
  }
}

void print_data(FILE *output, struct binstats *bs, int64_t count) {
  fprintf(output, "#Distance MedianDV D+ D-\n");
  fprintf(output, "#For halos between %g-%g Msun (no h)\n", MIN_MASS, MAX_MASS);
  fprintf(output, "#Num halos: %"PRId64"\n", count);
  print_medians(bs, output);
}

int main(int argc, char **argv) {
  int64_t i, j, p_next, count=0, p_start;
  struct bgc2_header hdr;
  struct binstats *bs_self, *bs_all;
  float bulkvel[3];
  char buffer[1024];
  FILE *output;
  char dirname[100];

  bs_self = init_binstats(BIN_START, BIN_END, BPDEX, 1, 1);
  bs_all = init_binstats(BIN_START, BIN_END, BPDEX, 1, 1);
  if (argc < 2) {
    printf("Usage: %s min_mass max_mass file1.bgc2 ...\n", argv[0]);
    exit(1);
  }
  MIN_MASS = atof(argv[1]);
  MAX_MASS = atof(argv[2]);
  assert(MIN_MASS>0 && MAX_MASS>MIN_MASS);

  for (i=3; i<argc; i++) {
    num_g = num_p = 0;
    load_bgc2(argv[i], &hdr, &grps, &num_g, &parts, &num_p);
    if (i==3) {
      snprintf(dirname, 100, "veldif_z%.2f_m%.2f_m%.2f", hdr.redshift, 
	       log10(MIN_MASS), log10(MAX_MASS));
      mkdir(dirname, 0755);
    }
    p_next = 0;
    for (j=0; j<num_g; j++) {
      p_next += grps[j].npart;
      if (grps[j].mass < MIN_MASS*hdr.Hubble0 || 
	  grps[j].mass > MAX_MASS*hdr.Hubble0 ||
	  grps[j].parent_id >= 0) continue;
      count++;
      p_start = p_next - grps[j].npart;
      calc_bulkvel(bulkvel, parts + p_start, grps[j].npart);
      sprintf(buffer, "%s/veldif_self_%05"PRId64".dat", dirname, count);
      output = check_fopen(buffer, "w");
      bin_velocities(parts + p_start, grps[j].npart_self, grps[j].pos,
		     bulkvel, bs_self, output);
      fclose(output);
      sprintf(buffer, "%s/veldif_all_%05"PRId64".dat", dirname, count);
      output = check_fopen(buffer, "w");
      bin_velocities(parts + p_start, grps[j].npart, grps[j].pos, bulkvel,
		     bs_all, output);
      fclose(output);
    }
  }

  sprintf(buffer, "%s/veldif_self.dat", dirname);
  output = check_fopen(buffer, "w");
  print_data(output, bs_self, count);
  fclose(output);
  sprintf(buffer, "%s/veldif_all.dat", dirname);
  output = check_fopen(buffer, "w");
  print_data(output, bs_all, count);
  fclose(output);
  return 0;
}



void load_bgc2(char *filename, struct bgc2_header *hdr,
	       GROUP_DATA_RMPVMAX **groups, int64_t *num_groups,
	       PARTICLE_DATA_PV **pdata, int64_t *num_parts)
{
  FILE *input;
  int64_t new_group_size, new_part_size;
  int64_t i, p_start;

  assert(sizeof(struct bgc2_header) == BGC2_HEADER_SIZE);
  input = check_fopen(filename, "rb");

  fread_fortran(hdr, BGC2_HEADER_SIZE, 1, input, 0);
  assert(hdr->magic == BGC_MAGIC);
  assert(hdr->version == 2);
  assert(hdr->format_group_data == GDATA_FORMAT_RMPVMAX);

  new_group_size = sizeof(GROUP_DATA_RMPVMAX)*((*num_groups)+hdr->ngroups);
  *groups = check_realloc(*groups, new_group_size, "Allocating groups.");
  fread_fortran((*groups) + (*num_groups), sizeof(GROUP_DATA_RMPVMAX), 
		hdr->ngroups, input, 0);
  *num_groups += hdr->ngroups;
  
  new_part_size = sizeof(PARTICLE_DATA_PV)*((*num_parts)+hdr->npart);
  *pdata = check_realloc(*pdata, new_part_size, "Allocating particles");
  p_start = 0;
  for (i=0; i<hdr->ngroups; i++) {
    fread_fortran((*pdata) + p_start, sizeof(PARTICLE_DATA_PV), 
		  groups[0][i].npart, input, 0);
    p_start += groups[0][i].npart;
  }
  *num_parts += hdr->npart;
  fclose(input);
}



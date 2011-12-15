#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "meta_io.h"
#include "../config_vars.h"
#include "../rockstar.h"
#include "../groupies.h"
#include "../check_syscalls.h"

void *output_buffer = NULL;
int64_t buffered = 0;

void fill_binary_header(struct binary_output_header *bh,
			int64_t snap, int64_t chunk) {
  bh->magic = ROCKSTAR_MAGIC;
  bh->snap = snap;
  bh->chunk = chunk;
  bh->scale = SCALE_NOW;
  bh->Om = Om;
  bh->Ol = Ol;
  bh->h0 = h0;
  bh->box_size = BOX_SIZE;
  bh->particle_mass = PARTICLE_MASS;  
}


void read_binary_header_config(struct binary_output_header *bh) {
  SCALE_NOW = bh->scale;
  Om = bh->Om;
  Ol = bh->Ol;
  h0 = bh->h0;
  BOX_SIZE = bh->box_size;
  PARTICLE_MASS = bh->particle_mass;
}

void output_particles_internal(int64_t snap, int64_t chunk) {
  char buffer[1024];
  struct binary_output_header bh;
  FILE *output;
  
  memset(&bh, 0, sizeof(struct binary_output_header));
  fill_binary_header(&bh, snap, chunk);
  bh.num_particles = num_p;
  bh.particle_type = PARTICLE_TYPE_FULL;

  get_output_filename(buffer, 1024, snap, chunk, "rbin");
  output = check_fopen(buffer, "wb");
  check_fwrite(&bh, sizeof(struct binary_output_header), 1, output);
  check_fwrite(p, sizeof(struct particle), num_p, output);
  fclose(output);
}

void load_particles_internal(char *filename, struct particle **part, int64_t *num_part) {
  FILE *input;
  struct binary_output_header bh;
  
  input = check_fopen(filename, "rb");
  check_fread(&bh, sizeof(struct binary_output_header), 1, input);
  assert(bh.magic == ROCKSTAR_MAGIC);
  assert(bh.particle_type == PARTICLE_TYPE_FULL);
  assert(bh.num_particles >= 0);
  read_binary_header_config(&bh);
  *part = check_realloc(*part, sizeof(struct particle)*bh.num_particles,
			"Allocating particles");
  check_fread(*part, sizeof(struct particle), bh.num_particles, input);
  *num_part = bh.num_particles;
  fclose(input);
}


inline void _clear_buffer(FILE *output) {
  check_fwrite(output_buffer, buffered, 1, output);
  buffered = 0;
}

inline void _append_to_buffer(void *src, int64_t size, FILE *output) {
  if (buffered + size > OUTPUT_BUFFER_SIZE) _clear_buffer(output);
  if (size > OUTPUT_BUFFER_SIZE) { check_fwrite(src, size, 1, output); return; }
  memcpy(output_buffer+buffered, src, size);
  buffered+=size;
}

void output_binary(int64_t id_offset, int64_t snap, int64_t chunk, float *bounds) {
  float max[3]={0}, min[3]={0};
  struct halo tmp;
  char buffer[1024];
  int64_t i,j, id=0;
  FILE *output;
  struct binary_output_header bheader;

  if (!output_buffer)
    output_buffer = check_realloc(output_buffer, OUTPUT_BUFFER_SIZE,
				  "Allocating output buffer");

  get_output_filename(buffer, 1024, snap, chunk, "bin");
  output = check_fopen(buffer, "wb");

  memset(&bheader, 0, sizeof(struct binary_output_header));
  _append_to_buffer(&bheader, sizeof(struct binary_output_header), output);

  //Output Halos
  if (num_halos) { 
    memcpy(min, halos[0].pos, sizeof(float)*3);
    memcpy(max, halos[0].pos, sizeof(float)*3);
  }
  for (i=0; i<num_halos; i++) {
    if (!_should_print(halos+i, bounds)) continue;
    for (j=0; j<3; j++) {
      if (min[j]>halos[i].pos[j]) min[j] = halos[i].pos[j];
      if (max[j]<halos[i].pos[j]) max[j] = halos[i].pos[j];
    }
    tmp = halos[i];
    tmp.id = id+id_offset;
    tmp.p_start = bheader.num_particles;
    _append_to_buffer(&tmp, sizeof(struct halo), output);
    bheader.num_halos++;
    bheader.num_particles+=tmp.num_p;
    id++;
  }

  //Output Particles
  for (i=0; i<num_halos; i++) {
    if (!_should_print(halos+i, bounds)) continue;
    for (j=0; j<halos[i].num_p; j++)
      _append_to_buffer(&(p[halos[i].p_start+j].id), sizeof(int64_t), output);
  }
  _clear_buffer(output);

  //Output header
  fill_binary_header(&bheader, snap, chunk);
  bheader.particle_type = PARTICLE_TYPE_IDS;
  if (bounds) memcpy(bheader.bounds, bounds, sizeof(float)*6);
  else { memcpy(bheader.bounds, min, sizeof(float)*3);
    memcpy(&(bheader.bounds[3]), max, sizeof(float)*3);
  }
  rewind(output);
  check_fwrite(&bheader, sizeof(struct binary_output_header), 1, output);
  fclose(output);  
}


void load_binary_header(int64_t snap, int64_t chunk, 
			struct binary_output_header *bheader) {
  char buffer[1024];
  FILE *input;
  get_output_filename(buffer, 1024, snap, chunk, "bin");
  input = check_fopen(buffer, "rb");
  check_fread(bheader, sizeof(struct binary_output_header), 1, input);
  assert(bheader->magic == ROCKSTAR_MAGIC);
  fclose(input);
}


void load_binary_halos(int64_t snap, int64_t chunk, 
      struct binary_output_header *bheader, struct halo **halos,
		       int64_t **part_ids)
{
  char buffer[1024];
  FILE *input;
  int64_t i, j;

  get_output_filename(buffer, 1024, snap, chunk, "bin");
  input = check_fopen(buffer, "rb");
  check_fread(bheader, sizeof(struct binary_output_header), 1, input);
  assert(bheader->magic == ROCKSTAR_MAGIC);
  assert(bheader->num_halos >= 0);
  assert(bheader->num_particles >= 0);
  *halos = check_realloc(*halos, sizeof(struct halo)*bheader->num_halos,
			 "Allocating room for halos.");
  *part_ids = check_realloc(*part_ids, sizeof(int64_t)*bheader->num_particles,
			    "Allocating room for particle IDs.");
  i = check_fread(*halos, sizeof(struct halo), bheader->num_halos, input);
  j = check_fread(*part_ids, sizeof(int64_t), bheader->num_particles, input);
  if (i!=bheader->num_halos || j!=bheader->num_particles) {
    fprintf(stderr, "[Error] Truncated input file %s!\n", buffer);
    exit(1);
  }
  fclose(input);

  read_binary_header_config(bheader);
}

#ifndef RESORT_OUTPUTS_H
#define RESORT_OUTPUTS_H

#include "gravitational_consistency.h"

void fork_and_print_halos(char *filename, struct halo_stash *h);
inline int64_t id_to_index(struct halo_stash h, int64_t id);
int sort_dead_by_mvir(const void *a, const void *b);
int sort_by_mvir(const void *a, const void *b);
int sort_by_desc(const void *a, const void *b);
void clear_halo_stash(struct halo_stash *h);
void zero_halo_stash(struct halo_stash *h);
void resort_halos(int last_output);

#endif /* RESORT_OUTPUTS_H */

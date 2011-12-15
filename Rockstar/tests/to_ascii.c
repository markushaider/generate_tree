#include <string.h>
#include "io/meta_io.h"
#include "groupies.h"

int main(int argc, char **argv) {
  struct binary_output_header bheader;
  int64_t *pids = NULL;
  load_binary_halos(0,0,&bheader, &halos, &pids);
  num_halos = bheader.num_halos;
  output_and_free_halos(0, 1, 0, NULL); 
  return 0;
}

#ifndef _CLIENT_H_
#define _CLIENT_H_

#define PARTICLE_RECIPIENT 0
#define HALO_RECIPIENT 1
#define PARTICLE_REALLOC_NUM 100000


struct recipient {
  int c;
  char *address;
  char *port;
  float bounds[6];
  void *buffer;
  int64_t buffered;
};

void client(void);
void recv_config(int c);
void send_config(int c);

#endif /* _CLIENT_H_ */

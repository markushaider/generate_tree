#ifndef _RSOCKET_H_
#define _RSOCKET_H_

#include <inttypes.h>
#define RECEIVER_FLAG 1
#define NEW_FLAG 2

struct rsocket {
  int fd;
  int64_t magic;
  int64_t flags;
  char *address;
  int port;
};

int start_server(char *host, int port, int fail);
struct rsocket refresh_rsocket(struct rsocket r);
struct rsocket connect_rsocket(char *address, int port);
struct rsocket reconnect_rsocket(struct rsocket r);
void close_rsocket(struct rsocket r);
int64_t send_to_rsocket(struct rsocket *r, void *data, int64_t length);
int64_t recv_from_rsocket(struct rsocket *r, void *data, int64_t length);
void *recv_and_alloc_r(struct rsocket *r, void *data, int64_t length);
int64_t send_msg_r(struct rsocket *r, void *data, int64_t length);
void *recv_msg_r(struct rsocket *r, void *data, int64_t *length, int64_t offset);
void *recv_msg_nolength_r(struct rsocket *r, void *data);


#endif /* _RSOCKET_H_ */

#ifndef _INET_SOCKET_H_
#define _INET_SOCKET_H_
#include <stdint.h>

#define RETRIES 10
#define TIMEOUT 2

int select_fd(int fd, int write, int timeout_secs);
int connect_to_addr(char *host, char *port);
int listen_at_addr(char *host, char *port);
int accept_connection(int s, char **address, int *port);
int64_t send_to_socket(int s, void *data, int64_t length);
int64_t recv_from_socket(int s, void *data, int64_t length);
void *recv_and_alloc(int s, void *data, int64_t length);
int64_t send_msg(int s, void *data, int64_t length);
void *recv_msg(int s, void *data, int64_t *length, int64_t offset);
void *recv_msg_nolength(int s, void *data);
void set_network_io_error_cb(void (*cb)(int), int data);
void *socket_check_realloc(void *ptr, size_t size, char *reason);
int get_recvbuf_size(int s);

#endif /* _INET_SOCKET_H_ */

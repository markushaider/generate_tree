#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include <assert.h>
#include "socket.h"
#include "rsocket.h"

#define LOCK pthread_mutex_lock(&socket_lock)
#define UNLOCK pthread_mutex_unlock(&socket_lock)

static struct rsocket *sockets = NULL;
static int64_t num_rsockets = 0;
static int64_t new_rsockets = 0;
static int server_port = 0;
static int hard_fail = 0;
static char *hostname = NULL;
pthread_t server_thread;
static pthread_mutex_t socket_lock = PTHREAD_MUTEX_INITIALIZER;

void ignore_error_cb(int ignore) { }

int64_t _find_socket_from_magic(uint64_t magic, int lock) {
  int64_t i;
  if (lock) LOCK;
  for (i=0; i<num_rsockets; i++)
    if (sockets[i].magic == magic) break;
  if (lock) UNLOCK;
  return i;
}

void *server_pthread(void *args) {
  int64_t i;
  int portnum, s=-1;
  char port[20];
  struct rsocket r;
  LOCK;
  for (i=0; i<5000&&s<0; i+=29) {
    portnum = server_port + i;
    snprintf(port, 10, "%d", portnum);
    s = listen_at_addr(hostname, port);
  }
  
  if (i>=5000) {
    if (hard_fail>0) {
      fprintf(stderr, "[Error] Couldn't start server at %s:%d-%d!\n",
	      hostname, server_port, portnum);
      exit(1);
    } else {
      hard_fail = -1;
      UNLOCK;
      pthread_exit(NULL);
    }
  }
  server_port = portnum;
  UNLOCK;

  r.address = 0;
  r.flags = RECEIVER_FLAG | NEW_FLAG;  
  while (1) {
    r.fd = accept_connection(s, &r.address, &r.port);
    recv_from_socket(r.fd, &r.magic, sizeof(uint64_t));
    LOCK;
    i = _find_socket_from_magic(r.magic, 0);
    if (i<num_rsockets) {
      close(sockets[i].fd);
      free(sockets[i].address);
      sockets[i] = r;
      sockets[i].flags -= sockets[i].flags & NEW_FLAG;
    } else {
      if (!(num_rsockets%10))
	sockets = socket_check_realloc(sockets, (num_rsockets+10)*
				       sizeof(struct rsocket), "Rsockets");
      sockets[num_rsockets] = r;
      num_rsockets++;
      new_rsockets++;
    }
    UNLOCK;
  }
  pthread_exit(NULL);
}


int start_server(char *host, int port, int fail) {
  LOCK;
  hostname = host;
  server_port = port;
  hard_fail = fail;
  UNLOCK;
  set_network_io_error_cb(ignore_error_cb, 0);
  return pthread_create(&server_thread, NULL, server_pthread, NULL);
}


struct rsocket refresh_rsocket(struct rsocket r)
{
  int64_t i = _find_socket_from_magic(r.magic, 1);
  if (i<num_rsockets) return sockets[i];
  return r;
}

struct rsocket connect_rsocket(char *address, int portnum) {
  struct rsocket r;
  char port[10];
  r.address = 0;
  r.port = portnum;
  snprintf(port, 10, "%d", r.port);
  r.flags = 0;
  r.magic = ((uint64_t)rand()) << (uint64_t)32;
  r.magic += rand();
  r.fd = connect_to_addr(address, port);
  if (r.fd < 0) return r;
  r.address = strdup(address);
  LOCK;
  if (!(num_rsockets%10))
    sockets = socket_check_realloc(sockets, (num_rsockets+10)*
				       sizeof(struct rsocket), "Rsockets");
  sockets[num_rsockets] = r;
  num_rsockets++;
  UNLOCK;
  return r;
}

struct rsocket reconnect_rsocket(struct rsocket r)
{
  char port[10];
  struct rsocket new_rsocket;
  if (r.flags & RECEIVER_FLAG) return refresh_rsocket(r);
  int64_t i = _find_socket_from_magic(r.magic, 1);

  snprintf(port, 10, "%d", r.port);
  if (i==num_rsockets) {
    new_rsocket = connect_rsocket(r.address, r.port);
    free(r.address);
    r = new_rsocket;
  } else {
    r.fd = connect_to_addr(r.address, port);
    LOCK;
    sockets[i].fd = r.fd;
    UNLOCK;
    if (r.fd > -1) send_to_socket(r.fd, &r.magic, sizeof(uint64_t));
  }
  return r;
}

void close_rsocket(struct rsocket r) {
  int64_t i = _find_socket_from_magic(r.magic, 1);
  if (i<num_rsockets) {
    LOCK;
    if (sockets[i].address) free(sockets[i].address);
    if (sockets[i].fd > -1) close(sockets[i].fd);
    sockets[i] = sockets[num_rsockets-1];
    num_rsockets--;
    UNLOCK;
  }
}


int64_t send_to_rsocket(struct rsocket *r, void *data, int64_t length) {
  int64_t sent;
  assert(r->fd >= 0); 
  while (r->fd < 0 || ((sent = send_to_socket(r->fd, data, length)) < 0)) {
    select_fd(-1, 0, 0);
    *r = reconnect_rsocket(*r);
  }
  return sent;
}

int64_t recv_from_rsocket(struct rsocket *r, void *data, int64_t length) {
  int64_t recvd;
  assert(r->fd >= 0); 
  while (r->fd < 0 || ((recvd = recv_from_socket(r->fd, data, length)) < 0)) {
    select_fd(-1, 0, 0);
    *r = reconnect_rsocket(*r);
  }
  return recvd;
}

void *recv_and_alloc_r(struct rsocket *r, void *data, int64_t length) {
  data = socket_check_realloc(data, length, "receive buffer");
  assert(data != NULL || length <= 0);
  recv_from_rsocket(r, data, length);
  return data;
}

int64_t send_msg_r(struct rsocket *r, void *data, int64_t length) {
  send_to_rsocket(r, &length, sizeof(int64_t));
  return send_to_rsocket(r, data, length);
}

void *recv_msg_r(struct rsocket *r, void *data, int64_t *length, int64_t offset) {
  int64_t incoming;
  if (recv_from_rsocket(r, &incoming, sizeof(int64_t))<0) return data;
  if (*length < (offset + incoming)) {
    data = socket_check_realloc(data, offset + incoming, "receive buffer");
    assert(data != NULL);
    *length = offset+incoming;
  }
  recv_from_rsocket(r, data+offset, incoming);
  return data;
}

void *recv_msg_nolength_r(struct rsocket *r, void *data) {
  int64_t length = 0;
  return(recv_msg_r(r, data, &length, 0));
}

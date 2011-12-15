#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/select.h>
#include <unistd.h>
#include <netinet/in.h>
#include <errno.h>
#include <netdb.h>
#include <assert.h>
#include <arpa/inet.h>
#include "socket.h"

void *socket_check_realloc(void *ptr, size_t size, char *reason) {
  void *res = realloc(ptr, size);
  if ((res == NULL) && (size > 0)) {
    fprintf(stderr, "[Error] Failed to allocate memory (%s)!\n", reason);
    exit(1);
  }
  return res;
}

//Using tutorials from http://beej.us/guide/bgnet/output/html/multipage/syscalls.html

void default_cb(int n) {
  perror("[Error] Network IO Failure");
  exit(1);
  //assert(0);
}

void (*io_error_cb)(int) = default_cb;
int io_error_data = 0;

void set_network_io_error_cb(void (*cb)(int), int data) {
  io_error_cb = cb;
  io_error_data = data;
}

struct addrinfo *default_addrinfo(char *host, char *port) {
  struct addrinfo hints = {0};
  struct addrinfo *res = NULL;
  int status;
  
  hints.ai_socktype = SOCK_STREAM;
  hints.ai_family = AF_UNSPEC;
  hints.ai_flags = (host) ? 0 : AI_PASSIVE;

  if (!(status = getaddrinfo(host, port, &hints, &res))) return res;
  fprintf(stderr, "[Error] Couldn't open %s:%s!  (Err: %s)\n", 
	  host, port, gai_strerror(status));
  exit(1);
}

int select_fd(int fd, int write, int timeout_secs) {
  fd_set fset;
  struct timeval timeout;
  FD_ZERO(&fset);
  timeout.tv_sec = timeout_secs;
  timeout.tv_usec = rand()%1000000;
  if (fd>-1) FD_SET(fd, &fset);
  if (write) return select(fd+1, NULL, &fset, NULL, &timeout);
  return select(fd+1, NULL, &fset, NULL, &timeout);
}

int connect_to_addr(char *host, char *port) {
  int s, i;
  struct addrinfo *res = default_addrinfo(host, port);
  s = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
  if (s < 0) {
    fprintf(stderr, "[Error] Couldn't open socket for address %s:%s!  (Err: %s)\n",
	    host, port, strerror(errno));
    exit(1);
  }

  for (i=0; i<RETRIES; i++) {
    if (!connect(s, res->ai_addr, res->ai_addrlen)) break;
    select_fd(-1, 0, 0);
  }
  if (i==RETRIES) {
    fprintf(stderr, "[Error] Failed to connect to %s:%s! (Err: %s;\n",
    	    host, port, strerror(errno));
    fprintf(stderr, " This error may mean that the connection was refused.)\n");
    exit(1);
  }
  freeaddrinfo(res);
  return s;
}

int listen_at_addr(char *host, char *port) {
  int s, i, set=1;
  struct addrinfo *res = default_addrinfo(host, port);
  s = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
  if (s < 0) {
    /*fprintf(stderr, "Couldn't open socket for address %s:%s!  (Err: %s)\n",
	    host, port, strerror(errno));
    exit(1);*/
    return -1;
  }

  if (setsockopt(s, SOL_SOCKET, SO_REUSEADDR, &set, sizeof(int))<0) {
    /*fprintf(stderr, "Couldn't set socket on %s:%s to reuse port!  (Err: %s)\n",
      host, port, strerror(errno));
      exit(1); */
    close(s);
    return -1;
  }

  for (i=0; i<RETRIES; i++)
    if (!bind(s, res->ai_addr, res->ai_addrlen)) break;

  if (i<RETRIES) 
    for (i=0; i<RETRIES; i++) 
      if (!listen(s, SOMAXCONN)) break;
  if (i==RETRIES) {
    /*fprintf(stderr, "Failed to accept connections at %s:%s! (Err: %s)\n",
	    host, port, strerror(errno));
    exit(1);*/
    close(s);
    return -1;
  }
  freeaddrinfo(res);
  return s;
}

int accept_connection(int s, char **address, int *port) {
  struct sockaddr_storage peer;
  socklen_t length = sizeof(peer);
  int a, size = (INET6_ADDRSTRLEN+1)*sizeof(char);
  a = accept(s, (struct sockaddr *)(&peer), &length);
  if (a<0) {
    perror("[Error] Connection accept failed");
    assert(0);
  }

  if (address) {
    if (*address) free(*address);
    *address = (char *)malloc(size);
  }

  if (peer.ss_family == AF_INET) {
    struct sockaddr_in *s = (struct sockaddr_in *)&peer;
    if (port) *port = ntohs(s->sin_port);
    if (address && *address) inet_ntop(AF_INET, &s->sin_addr, *address, size);
  } else { // Assume AF_INET6
    struct sockaddr_in6 *s = (struct sockaddr_in6 *)&peer;
    if (port) *port = ntohs(s->sin6_port);
    if (address && *address) inet_ntop(AF_INET6, &s->sin6_addr, *address, size);
  }
  return a;
}

int64_t send_to_socket(int s, void *data, int64_t length) {
  int64_t sent = 0, n = 0;
  while (sent < length) {
    n = send(s, data+sent, length-sent, 0);
    if (n<0) {
      io_error_cb(io_error_data);
      return -1;
    }
    sent+=n;
  }
  return (sent);
}

int64_t recv_from_socket(int s, void *data, int64_t length) {
  int64_t received = 0, n = 0;
  while (received < length) {
    n = recv(s, data+received, length-received, 0);
    if (n<0) {
      io_error_cb(io_error_data);
      return -1;
    }
    received+=n;
  }
  return (received);
}

void *recv_and_alloc(int s, void *data, int64_t length) {
  data = socket_check_realloc(data, length, "receive buffer");
  assert(data != NULL || length <= 0);
  recv_from_socket(s, data, length);
  return data;
}

int64_t send_msg(int s, void *data, int64_t length) {
  send_to_socket(s, &length, sizeof(int64_t));
  return send_to_socket(s, data, length);
}

void *recv_msg(int s, void *data, int64_t *length, int64_t offset) {
  int64_t incoming;
  if (recv_from_socket(s, &incoming, sizeof(int64_t))<0) return data;
  if (*length < (offset + incoming)) {
    data = socket_check_realloc(data, offset + incoming, "receive buffer");
    assert(data != NULL);
    *length = offset+incoming;
  }
  recv_from_socket(s, data+offset, incoming);
  return data;
}

void *recv_msg_nolength(int s, void *data) {
  int64_t length = 0;
  return(recv_msg(s, data, &length, 0));
}

int get_recvbuf_size(int s) {
  socklen_t optlen = sizeof(int);
  int recvbuf;
  getsockopt(s, SOL_SOCKET, SO_RCVBUF, &recvbuf, &optlen);
  return recvbuf;
}

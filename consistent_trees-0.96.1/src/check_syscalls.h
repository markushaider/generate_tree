#ifndef CHECK_SYSCALLS_H
#define CHECK_SYSCALLS_H

FILE *check_fopen(char *filename, char *mode);
void *check_realloc(void *ptr, size_t size, char *reason);

#endif /* CHECK_SYSCALLS_H */

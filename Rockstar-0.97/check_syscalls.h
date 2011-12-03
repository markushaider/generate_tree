#ifndef CHECK_SYSCALLS_H
#define CHECK_SYSCALLS_H

void system_error(char *errmsg);
FILE *check_fopen(char *filename, char *mode);
void *check_realloc(void *ptr, size_t size, char *reason);
size_t check_fread(void *ptr, size_t size, size_t nitems, FILE *stream);
size_t check_fwrite(void *ptr, size_t size, size_t nitems, FILE *stream);
void check_fseeko(FILE *stream, off_t offset, int whence);

#endif /* CHECK_SYSCALLS_H */

CFLAGS=-m64 -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -lm -Wall
OFLAGS = -O3 -std=c99
DEBUGFLAGS = -g -O -std=c99
PROFFLAGS = -g -pg -O2 -std=c99
CC = gcc

all:
	@make reg EXTRA_FLAGS="$(OFLAGS)"

debug:
	@make reg EXTRA_FLAGS="$(DEBUGFLAGS)"

prof:
	@make reg EXTRA_FLAGS="$(PROFFLAGS)"


reg:
	$(CC) $(CFLAGS) example.c read_tree.c $(EXTRA_FLAGS) stringparse.c check_syscalls.c -o read_tree

clean:
	rm -f *~
	rm -f read_tree
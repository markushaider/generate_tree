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
	$(CC) last_minor_m.c ../read_tree/stringparse.c ../read_tree/check_syscalls.c ../src/universe_time.c $(CFLAGS) $(EXTRA_FLAGS) -o last_minor_m
	$(CC) find_parents.c stringparse.c ../src/check_syscalls.c  -o find_parents $(CFLAGS) $(EXTRA_FLAGS)
	$(CC) subhalo_stats_orig.c stringparse.c ../src/check_syscalls.c  ../src/litehash.c -o subhalo_stats_orig $(CFLAGS) $(EXTRA_FLAGS)
	$(CC) subhalo_destroyed.c stringparse.c ../src/check_syscalls.c  ../src/litehash.c -o subhalo_destroyed $(CFLAGS) $(EXTRA_FLAGS)
	$(CC) subhalo_stats.c stringparse.c ../src/check_syscalls.c  ../src/litehash.c -o subhalo_stats $(CFLAGS) $(EXTRA_FLAGS)
	$(CC) mass_function.c stringparse.c ../src/halo_io.c ../src/check_syscalls.c ../src/masses.c $(CFLAGS) $(EXTRA_FLAGS) -o mass_function
	$(CC) subhalo_mistracked.c ../read_tree/read_tree.c ../read_tree/stringparse.c ../read_tree/check_syscalls.c $(CFLAGS) $(EXTRA_FLAGS) -o subhalo_mistracked
	$(CC) rpeak.c ../read_tree/stringparse.c ../read_tree/check_syscalls.c ../src/universe_time.c $(CFLAGS) $(EXTRA_FLAGS) -o rpeak
	$(CC) sub_merged.c ../read_tree/stringparse.c ../read_tree/check_syscalls.c ../src/universe_time.c $(CFLAGS) $(EXTRA_FLAGS) -o sub_merged

CFLAGS=-m64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_BSD_SOURCE -D_POSIX_SOURCE -D_POSIX_C_SOURCE=200809L -D_SVID_SOURCE -D_DARWIN_C_SOURCE -lm -Wall
SUPPORT=masses.c distance.c universe_time.c stringparse.c halo_io.c check_syscalls.c
EXTRA_SUPPORT=$(SUPPORT) read_config.c grav_config.c
OFLAGS = -O3 -std=c99 -lm
DEBUGFLAGS = -g -O -std=c99
PROFFLAGS = -g -pg -O2 -std=c99
CC = gcc

all:
	@make reg EXTRA_FLAGS="$(OFLAGS)"

no-periodic:
	@make reg EXTRA_FLAGS="$(OFLAGS) -DNO_PERIODIC" NP="_no_periodic"

debug:
	@make reg EXTRA_FLAGS="$(DEBUGFLAGS)"

prof:
	@make reg EXTRA_FLAGS="$(PROFFLAGS)"


reg:
	$(CC) $(CFLAGS) halo_evolve.c halo_evolve_lib.c $(SUPPORT) -o ../halo_evolve$(NP)  $(EXTRA_FLAGS)
	$(CC) $(CFLAGS) tidal.c tidal_lib.c $(EXTRA_SUPPORT) -o ../tidal$(NP)  $(EXTRA_FLAGS)
	$(CC) $(CFLAGS) gravitational_consistency.c gravitational_statistics.c halo_evolve_lib.c tidal_lib.c $(EXTRA_SUPPORT) -o ../gravitational_consistency$(NP) $(EXTRA_FLAGS)
	$(CC) $(CFLAGS) find_parents_and_cleanup.c cleanup_statistics.c litehash.c $(EXTRA_SUPPORT) -o ../find_parents_and_cleanup$(NP) $(EXTRA_FLAGS)
	$(CC) $(CFLAGS) resort_outputs.c $(EXTRA_SUPPORT) -o ../resort_outputs $(EXTRA_FLAGS) 
	$(CC) $(CFLAGS) assemble_halo_trees.c cached_io.c litehash.c $(EXTRA_SUPPORT) -o ../assemble_halo_trees $(EXTRA_FLAGS)

clean:
	rm -f *~
	rm -f ../find_parents_and_cleanup ../gravitational_consistency ../tidal ../halo_evolve ../resort_outputs ../halo_evolve_no_periodic ../gravitational_consistency_no_periodic ../find_parents_and_cleanup_no_periodic ../assemble_halo_trees
	rm -rf ../find_parents_and_cleanup.* ../gravitational_consistency.* ../tidal.* ../halo_evolve.* ../resort_outputs.*

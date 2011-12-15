//string(a,b) real(a,b) integer(a,b)
//Sets strings, floats, and integers; 
//     "a" is the variable name; 
//     "b" is the default.

string(FILE_FORMAT, "GADGET2");
real(PARTICLE_MASS, 0); //Auto-generated or auto-read

string(MASS_DEFINITION, "vir");
integer(MIN_HALO_OUTPUT_SIZE, 20);
real(FORCE_RES, 0.003); //In Mpc/h

real(SCALE_NOW, 1.0);
real(h0, 0.7);
real(Ol, 0.73);
real(Om, 0.27);

integer(GADGET_ID_BYTES, 4);
real(GADGET_MASS_CONVERSION, 1e10);
real(GADGET_LENGTH_CONVERSION, 1.0);
integer(GADGET_SKIP_NON_HALO_PARTICLES, 1);
integer(RESCALE_PARTICLE_MASS, 0);

real(TIPSY_LENGTH_CONVERSION, 1.0);
real(TIPSY_VELOCITY_CONVERSION, 1.0);

integer(PARALLEL_IO, 0);
string(PARALLEL_IO_SERVER_ADDRESS, "auto");
string(PARALLEL_IO_SERVER_PORT, "auto");
integer(PARALLEL_IO_WRITER_PORT, 32001);
string(PARALLEL_IO_SERVER_INTERFACE, "");
string(RUN_ON_SUCCESS, "");

string(INBASE, ".");
string(FILENAME,"tests/halo_nfw");
integer(STARTING_SNAP, 0);
integer(NUM_SNAPS, 1);
integer(NUM_BLOCKS, 1);
integer(NUM_READERS, 0);
integer(PRELOAD_PARTICLES, 0);
string(SNAPSHOT_NAMES, "");
string(LIGHTCONE_ALT_SNAPS, "");
string(BLOCK_NAMES, "");

string(OUTBASE, ".");
real(OVERLAP_LENGTH, 3.0);
integer(NUM_WRITERS, 1);
integer(FORK_READERS_FROM_WRITERS, 0);
integer(FORK_PROCESSORS_PER_MACHINE, 1);

string(OUTPUT_FORMAT, "BOTH");
integer(DELETE_BINARY_OUTPUT_AFTER_FINISHED, 0);
integer(FULL_PARTICLE_CHUNKS, 0);
string(BGC2_SNAPNAMES, "");

integer(BOUND_PROPS, 1);
integer(BOUND_OUT_TO_HALO_EDGE, 0);
integer(DO_MERGER_TREE_ONLY, 0);
integer(IGNORE_PARTICLE_IDS, 0);
real(TRIM_OVERLAP, 0);
real(ROUND_AFTER_TRIM, 1);
integer(LIGHTCONE, 0);
integer(PERIODIC, 1);

real3(LIGHTCONE_ORIGIN, "0 0 0");
real3(LIGHTCONE_ALT_ORIGIN, "0 0 0");

real3(LIMIT_CENTER, "0 0 0");
real(LIMIT_RADIUS, 0);

integer(SWAP_ENDIANNESS, 0);
integer(GADGET_VARIANT, 0);

real(FOF_FRACTION, 0.7);
real(FOF_LINKING_LENGTH, 0.28);
real(INCLUDE_HOST_POTENTIAL_RATIO, 0.3);
real(DOUBLE_COUNT_SUBHALO_MASS_RATIO, 1.0);
integer(TEMPORAL_HALO_FINDING, 1);
integer(MIN_HALO_PARTICLES, 10);
real(UNBOUND_THRESHOLD, 0.5);
integer(ALT_NFW_METRIC, 0);

integer(TOTAL_PARTICLES, 8589934592);
real(BOX_SIZE, 250); //In Mpc/h
integer(OUTPUT_HMAD, 0);
integer(OUTPUT_PARTICLES, 0);
integer(OUTPUT_LEVELS, 0);
real3(DUMP_PARTICLES, "0 0 0");

real(AVG_PARTICLE_SPACING, 0); //Auto-generated
integer(SINGLE_SNAP, 0);

#ifndef _GRAVITATIONAL_CONSISTENCY_VARS_H_
#define _GRAVITATIONAL_CONSISTENCY_VARS_H_

float Om=0.27;
float h0=0.70;
float Ol=0.73;

char *SCALEFILE = "/Volumes/Peter 1/Bolshoi/DescScales.txt";
char *INBASE = "/Volumes/Peter 1/Bolshoi/Redshift_Lists";
char *OUTBASE = "/Volumes/Peter 1/Bolshoi/Redshift_Lists_Bak";
char *TREE_OUTBASE = "/Volumes/Peter 1/Bolshoi/Trees";

float MAJOR_MERGER=0.3;
float MIN_MMP_MASS_RATIO=0.5;
float MIN_MMP_VMAX_RATIO=0.7;

float BOX_WIDTH = 250;
float BOX_DIVISIONS = 5;

float SOFTENING_LENGTH = 1;

int PADDING_TIMESTEPS=1; //Don't kill halos if they have link problems during the last 1 timestep.

int MIN_TIMESTEPS_SUB_TRACKED = 20; // Only keep tracks if they have been around for longer than this time.
int MIN_TIMESTEPS_SUB_MMP_TRACKED = 10; // Only keep tracks if they have been around for longer than this time.
int MIN_TIMESTEPS_TRACKED = 5; // Only keep tracks if they have been around for longer than this time.
float MAX_PHANTOM_FRACTION = 0.2; // Discard tracks where there are more than 20% phantoms.

float LAST_DITCH_SEARCH_LIMIT=1.0; /* For connecting halos which have "moved" up to this amount
				       times their virial radius */
float LAST_DITCH_VMAX_RATIO_1=1.1; /* For connecting halos which have "moved" unphysical amounts */
float LAST_DITCH_VMAX_RATIO_2=2.5; /* For connecting halos which have "moved" unphysical amounts */

int MAX_PHANTOM=4; /* max timesteps to keep phantom halo */
int MAX_PHANTOM_SMALL=2; /* max timesteps to keep small phantom halo */
int SMALL_PARTICLE_LIMIT=49; /* Halos smaller than this size get kept around for less time. */
/* SHOULD NOT BE SET TO MORE THAN 50 for BOLSHOI and BDM results. */

float TIDAL_FORCE_LIMIT=0.1;
int RECURSION_LIMIT=5;
float METRIC_LIMIT=7;
float UNPHYSICAL=22; //Should be (3*METRIC_LIMIT+1);
float METRIC_BREAK_LIMIT=3.2; //Below which we break a link.
float MASS_RES_OK=1e11; //Halo mass above which there are probably not resolution issues.

#endif /* _GRAVITATIONAL_CONSISTENCY_VARS_H_ */

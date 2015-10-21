/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

*/
#ifndef __MATVEC_H
#include "matvec.h"
#endif

#ifndef __GENCLASSES_H
#include "genclasses.h"
#endif


#ifndef __DISTANCES_H
#define  __DISTANCES_H

/*
   ==================
   === STRUCTURES ===
   ==================
*/


typedef struct{
	mat_int *mutation1, *mutation2, *nbcommon;
	int n;
} dna_dist;


typedef struct{
	mat_double *dist;
	int n;
} spatial_dist;


/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

dna_dist * alloc_dna_dist(int n);

spatial_dist * alloc_spatial_dist(int n);



/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

void free_dna_dist(dna_dist * in);

void free_spatial_dist(spatial_dist * in);





/*
   =================
   === AUXILIARY ===
   =================
*/

bool is_atgc(char in);

int get_mutation1(dna_dist * in, int i, int j);

int get_mutation2(dna_dist * in, int i, int j);

int get_nbcommon(dna_dist * in, int i, int j);

double get_spatial_dist(spatial_dist * in, int i, int j);



/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/

void print_dna_dist(dna_dist *in);

void print_spatial_dist(spatial_dist *in);

dna_dist * compute_dna_distances(list_dnaseq *in, int mut_model);

spatial_dist * doublevec2spatial_dist(double *in, int *n);


#endif

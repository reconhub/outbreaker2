/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), January 2012.
  Licence: GPL >=2.

*/

#include "common.h"
#include "matvec.h"
#include "genclasses.h"
#include "distances.h"




/*
  ====================
  === CONSTRUCTORS ===
  ====================
*/

/* ALLOC A DNA_DIST OBJECT */
dna_dist * alloc_dna_dist(int n){
    dna_dist * out = (dna_dist *) malloc(sizeof(dna_dist));
    if(out==NULL){
      error("\n[in: distances.c->alloc_dna_dist]\nNo memory left for creating distance matrix. Exiting.\n");
	/* fprintf(stderr, "\n[in: distances.c->alloc_dna_dist]\nNo memory left for creating distance matrix. Exiting.\n"); */
	/* exit(1); */
    }

    out->mutation1 = alloc_mat_int(n,n);
    out->mutation2 = alloc_mat_int(n,n);
    out->nbcommon = alloc_mat_int(n,n);
    out->n = n;

    return out;
}



/* ALLOC A SPATIAL_DIST OBJECT */
spatial_dist * alloc_spatial_dist(int n){
    spatial_dist * out = (spatial_dist *) malloc(sizeof(spatial_dist));
    if(out==NULL){
      error("\n[in: distances.c->alloc_spatial_dist]\nNo memory left for creating distance matrix. Exiting.\n");
	/* fprintf(stderr, "\n[in: distances.c->alloc_spatial_dist]\nNo memory left for creating distance matrix. Exiting.\n"); */
	/* exit(1); */
    }

    out->dist = alloc_mat_double(n,n);
    out->n = n;

    return out;
}






/*
  ===================
  === DESTRUCTORS ===
  ===================
*/

void free_dna_dist(dna_dist * in){
    free_mat_int(in->mutation1);
    free_mat_int(in->mutation2);
    free_mat_int(in->nbcommon);
    free(in);
}



void free_spatial_dist(spatial_dist * in){
    free_mat_double(in->dist);
    free(in);
}





/*
  =================
  === AUXILIARY ===
  =================
*/

bool is_atgc(char in){
    if(in=='a' || in=='t' || in=='g' || in=='c') return TRUE;
    return FALSE;
}


int get_mutation1(dna_dist * in, int i, int j){
    return in->mutation1->rows[i]->values[j];
}


int get_mutation2(dna_dist * in, int i, int j){
    return in->mutation2->rows[i]->values[j];
}


int get_nbcommon(dna_dist * in, int i, int j){
    return in->nbcommon->rows[i]->values[j];
}


double get_spatial_dist(spatial_dist * in, int i, int j){
    return in->dist->rows[i]->values[j];
}


/*
  ===============================
  === MAIN EXTERNAL FUNCTIONS ===
  ===============================
*/

void print_dna_dist(dna_dist *in){
    Rprintf("\n - transitions -");
    print_mat_int(in->mutation1);

    Rprintf("\n - transversions -");
    print_mat_int(in->mutation2);

    Rprintf("\n - common nucleotides -");
    print_mat_int(in->nbcommon);
    Rprintf("\n");
}





void print_spatial_dist(spatial_dist *in){
    Rprintf("\n - distance -");
    print_mat_double(in->dist);
    Rprintf("\n");
}




/* GET ALL PAIRWISE DISTANCES (TRANSITIONS/TRANSVERSIONS) IN A LIST OF SEQUENCES */
dna_dist * compute_dna_distances(list_dnaseq *in, int mut_model){
    int i,j,k;
    int N=in->n, L=in->length;

    /* CREATE OUTPUT */
    dna_dist *out = alloc_dna_dist(N);

    /* CHECK MODEL */
    if(mut_model<0 || mut_model>2) error("\n[in: distances.c->compute_dna_distances]\nModel %d unknown.\n", mut_model);


    /* COMPUTE DISTANCES */
    /* MODEL 1: distance = nb of mutations */
    if(mut_model==1){
	/* for all unique pairs of sequences */
	for(i=0;i<(N-1);i++){
	    for(j=i+1;j<N;j++){
		/* for all pairs of nucleotides */
		for(k=0;k<L;k++){
		    if(is_atgc(in->list[i]->seq[k]) && is_atgc(in->list[j]->seq[k])){ /*if non-missing data*/
			/* one more nucleotide was comparable */
			out->nbcommon->rows[i]->values[j] = out->nbcommon->rows[i]->values[j] + 1;
			if(in->list[i]->seq[k] != in->list[j]->seq[k]){
				out->mutation1->rows[i]->values[j] = out->mutation1->rows[i]->values[j] + 1;
			}
		    } /* end if non-missing data*/
		} /* end for k */

		/* FILL IN THE SECOND HALF OF THE 'MATRIX' */
		out->mutation1->rows[j]->values[i] = out->mutation1->rows[i]->values[j];
		out->nbcommon->rows[j]->values[i] = out->nbcommon->rows[i]->values[j];
	    } /* end for j */
	} /* end for i */
    }

    /* MODEL 2: mutation1=transitions, mutation2=transversions */
    if(mut_model==2){
	/* for all unique pairs of sequences */
	for(i=0;i<(N-1);i++){
	    for(j=i+1;j<N;j++){
		/* for all pairs of nucleotides */
		for(k=0;k<L;k++){
		    if(is_atgc(in->list[i]->seq[k]) && is_atgc(in->list[j]->seq[k])){ /*if non-missing data*/
			/* one more nucleotide was comparable */
			out->nbcommon->rows[i]->values[j] = out->nbcommon->rows[i]->values[j] + 1;
			if(in->list[i]->seq[k] != in->list[j]->seq[k]){
			    /* transitions */
			    if((in->list[i]->seq[k]=='a' && in->list[j]->seq[k]=='g')
			       || (in->list[i]->seq[k]=='g' && in->list[j]->seq[k]=='a')
			       || (in->list[i]->seq[k]=='c' && in->list[j]->seq[k]=='t')
			       || (in->list[i]->seq[k]=='t' && in->list[j]->seq[k]=='c')) {
				out->mutation1->rows[i]->values[j] = out->mutation1->rows[i]->values[j] + 1;
			    } else { /* else it is a transversion*/
				out->mutation2->rows[i]->values[j] = out->mutation2->rows[i]->values[j] + 1;
			    }
			}
		    } /* end if non-missing data*/
		} /* end for k */

		/* FILL IN THE SECOND HALF OF THE 'MATRIX' */
		out->mutation1->rows[j]->values[i] = out->mutation1->rows[i]->values[j];
		out->mutation2->rows[j]->values[i] = out->mutation2->rows[i]->values[j];
		out->nbcommon->rows[j]->values[i] = out->nbcommon->rows[i]->values[j];
	    } /* end for j */
	} /* end for i */
    }

    /* SEQUENCES HAVE L NUCLEOTIDES IN COMMON WITH THEMSELVES */
    for(i=0;i<N;i++){
	out->nbcommon->rows[i]->values[i] = L;
    }

    /* RETURN RESULT */
    return out;
} /* end compute_dna_distances */





/* CREATE A SPATIAL_DIST OBJECT FROM AN INPUT FROM R */
spatial_dist * doublevec2spatial_dist(double *in, int *n){
    int i, j, count=0;

    /* ALLOC OUTPUT */
    spatial_dist *out = alloc_spatial_dist(*n);

  /* FILL IN THE OUTPUT */
    for(j=0;j<*n;j++){
	for(i=0;i<*n;i++){
	    out->dist->rows[i]->values[j] = in[count++];
	}
    }

    return out;
} /* end  doublevec2spatial_dist */


/*
  =========================
  === TESTING FUNCTIONS ===
  =========================
*/


/* int main(){ */
/* 	const int N=5, L=10; */
/* 	int i,j; */

/* 	/\* create a list of sequences *\/ */
/* 	list_dnaseq * test = alloc_list_dnaseq(N, L); */

/* 	for(i=0;i<N;i++){ */
/* 		for(j=0;j<L;j++){ */
/* 			if(i*j % 5 ==0) test->list[i]->seq[j] = 'a'; */
/* 			else if(i*j % 3 ==0)test->list[i]->seq[j] = 't'; */
/* 			else if(i*j % 2 ==0)test->list[i]->seq[j] = 'g'; */
/* 			else test->list[i]->seq[j] = 'c'; */
/* 		} */
/* 	} */

/* 	for(i=5;i<L;i++) */
/* 		test->list[0]->seq[i] = '-'; */
/* 	for(i=0;i<5;i++) */
/* 		test->list[N-1]->seq[i] = '-'; */

/* 	print_list_dnaseq(test); */

/* 	/\* compute distances *\/ */
/* 	dna_dist *out = compute_dna_distances(test); */


/* 	print_dna_dist(out); */

/* 	/\* free and return *\/ */
/* 	free_list_dnaseq(test); */
/* 	free_dna_dist(out); */

/* 	return 0; */
/* } */



/*
  gcc instructions

  gcc -o distances matvec.c genclasses.c distances.c -lgsl -lgslcblas -g

  ./distances

  valgrind --leak-check=full -v distances

*/

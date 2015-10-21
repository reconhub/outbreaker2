/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), January 2012.
  Licence: GPL >=2.

*/

#include "common.h"
#include "matvec.h"



/*
  ====================
  === CONSTRUCTORS ===
  ====================
*/

/* ALLOC A VECTOR OF INTEGERS OF SIZE N */
vec_int * alloc_vec_int(int n){
    vec_int *out = (vec_int *) malloc(sizeof(vec_int));
    if(out == NULL){
      error("\n[in: matvec.c->alloc_vec_int]\nNo memory left for creating vector of integers. Exiting.\n");
      /* fprintf(stderr, "\n[in: matvec.c->alloc_vec_int]\nNo memory left for creating vector of integers. Exiting.\n"); */
      /* exit(1); */
    }

    /* NOTE out->values is not allocated when n=0 */
    if(n>0){
	out->values = (int *) calloc(n, sizeof(int));
	if(out->values == NULL){
	  error("\n[in: matvec.c->alloc_vec_int]\nNo memory left for creating vector of integers. Exiting.\n");
	  /* fprintf(stderr, "\n[in: matvec.c->alloc_vec_int]\nNo memory left for creating vector of integers. Exiting.\n"); */
	  /* exit(1); */
	}
    }

    out->length = n;

    return(out);
}





/* ALLOC A VECTOR OF DOUBLEEGERS OF SIZE N */
vec_double * alloc_vec_double(int n){
    vec_double *out = (vec_double *) malloc(sizeof(vec_double));
    if(out == NULL){
      error("\n[in: matvec.c->alloc_vec_double]\nNo memory left for creating vector of double. Exiting.\n");
      /* fprintf(stderr, "\n[in: matvec.c->alloc_vec_double]\nNo memory left for creating vector of double. Exiting.\n"); */
      /* exit(1); */
    }

    /* NOTE out->values is not allocated when n=0 */
    if(n>0){
	out->values = (double *) calloc(n, sizeof(double));
	if(out->values == NULL){
	  error("\n[in: matvec.c->alloc_vec_double]\nNo memory left for creating vector of double. Exiting.\n");
	  /* fprintf(stderr, "\n[in: matvec.c->alloc_vec_double]\nNo memory left for creating vector of double. Exiting.\n"); */
	  /* exit(1); */
	}
    }

    out->length = n;

    return(out);
}





/* ALLOC EMPTY MAT_INT BETWEEN N OBJECTS */
/* (values initialized to 0) */
/* n rows, p columns */
mat_int * alloc_mat_int(int n, int p){
    int i;
    mat_int *out;

    /* allocate output */
    out = (mat_int *) malloc(sizeof(mat_int));
    if(out == NULL){
      error("\n[in: matvec.c->alloc_mat_int]\nNo memory left for creating distance matrix. Exiting.\n");
      /* fprintf(stderr, "\n[in: matvec.c->alloc_mat_int]\nNo memory left for creating distance matrix. Exiting.\n"); */
      /* exit(1); */
    }

    /* fill in content */
    out->rows = (vec_int **) calloc(n, sizeof(vec_int *));
    if(out->rows == NULL){
      error("\n[in: matvec.c->alloc_mat_int]\nNo memory left for creating distance matrix. Exiting.\n");
      /* fprintf(stderr, "\n[in: matvec.c->alloc_mat_int]\nNo memory left for creating distance matrix. Exiting.\n"); */
      /* exit(1); */
    }

    for(i=0;i<n;i++){
	out->rows[i] = alloc_vec_int(p);
    }

    out->n = n;
    out->p = p;

    /* return */
    return out;
}





/* ALLOC EMPTY MAT_DOUBLE BETWEEN N OBJECTS */
/* (values initialized to 0) */
mat_double * alloc_mat_double(int n, int p){
    int i;
    mat_double *out;

    /* allocate output */
    out = (mat_double *) malloc(sizeof(mat_double));
    if(out == NULL){
      error("\n[in: matvec.c->alloc_mat_double]\nNo memory left for creating distance matrix. Exiting.\n");
      /* fprintf(stderr, "\n[in: matvec.c->alloc_mat_double]\nNo memory left for creating distance matrix. Exiting.\n"); */
      /* exit(1); */
    }

    /* fill in content */
    out->rows = (vec_double **) calloc(n, sizeof(vec_double *));
    if(out->rows == NULL){
      error("\n[in: matvec.c->alloc_mat_double]\nNo memory left for creating distance matrix. Exiting.\n");
      /* fprintf(stderr, "\n[in: matvec.c->alloc_mat_double]\nNo memory left for creating distance matrix. Exiting.\n"); */
      /* exit(1); */
    }

    for(i=0;i<n;i++){
	out->rows[i] = alloc_vec_double(p);
    }

    out->n = n;
    out->p = p;

    /* return */
    return out;
}




/*
  ===================
  === DESTRUCTORS ===
  ===================
*/


void free_vec_int(vec_int *in){
    if(in->length > 0) free(in->values);
    free(in);
}


void free_mat_int(mat_int *in){
    int i;
    if(in->n > 0) {
	for(i=0;i<in->n;i++)
	    free_vec_int(in->rows[i]);
    }
    free(in->rows);
    free(in);
}




void free_vec_double(vec_double *in){
    if(in->length > 0) free(in->values);
    free(in);
}


void free_mat_double(mat_double *in){
    int i;
    if(in->n > 0) {
	for(i=0;i<in->n;i++)
	    free_vec_double(in->rows[i]);
    }
    free(in->rows);
    free(in);
}




/*
  ===============================
  === MAIN EXTERNAL FUNCTIONS ===
  ===============================
*/

int vec_int_i(vec_int *in, int i){
    if(i >= in->length || i<0) {
      error("\n[in: vec_int_i] Trying to access value %d in a vector of size %d\n",i,in->length);
      /* fprintf(stderr, "\nTrying to access value %d in a vector of size %d\n",i,in->length); */
      /* exit(1); */
    }
    return in->values[i];
}




int mat_int_ij(mat_int *in, int i, int j){
    if(i >= in->n || i<0 || j>=in->p || j<0) {
      error("\n[in: mat_int_i] Trying to access item (%d,%d) in a matrix of dimensions (%d,%d)\n", i, j, in->n, in->p);
      /* fprintf(stderr, "\nTrying to access item %d in a list of size %d\n",i,in->n); */
      /* exit(1); */
    }
    return vec_int_i(in->rows[i], j);
}





double vec_double_i(vec_double *in, int i){
    if(i >= in->length || i<0) {
      error("\n[in: vec_double_i] Trying to access value %d in a vector of size %d\n",i,in->length);
      /* fprintf(stderr, "\nTrying to access value %d in a vector of size %d\n",i,in->length); */
      /* exit(1); */
    }
    return in->values[i];
}




double mat_double_ij(mat_double *in, int i, int j){
  if(i >= in->n || i<0 || j>=in->p || j<0) {
    error("\n[in: mat_double_i] Trying to access item (%d,%d) in a matrix of dimensions (%d,%d)\n", i, j, in->n, in->p);
    /* fprintf(stderr, "\nTrying to access item %d in a list of size %d\n",i,in->n); */
    /* exit(1); */
  }
  return vec_double_i(in->rows[i], j);
}






/* print method */
void print_vec_int(vec_int *in){
    int i;
    Rprintf("Vector of %d integers: ", in->length);
    /* for(i=0;i<in->length;i++) printf("%d ", in->values[i]); */
    for(i=0;i<in->length;i++) Rprintf("%d ", vec_int_i(in,i));
    Rprintf("\n");
}




/* print method */
void print_mat_int(mat_int *in){
    int i,j;
    Rprintf("\n%dx%d matrix of integers",in->n,in->p);
    for(i=0;i<in->n;i++){
	Rprintf("\n");
	for(j=0;j<in->p;j++)
	    /* Rprintf("%d ", in->rows[i]->values[j]); */
	    Rprintf("%d ", mat_int_ij(in,i,j));
    }
    Rprintf("\n");
}






/* print method */
void print_vec_double(vec_double *in){
    int i;
    Rprintf("Vector of %d doubles: ", in->length);
    /* for(i=0;i<in->length;i++) printf("%d ", in->values[i]); */
    for(i=0;i<in->length;i++) Rprintf("%.3f ", vec_double_i(in,i));
    Rprintf("\n");
}





/* print method */
void print_mat_double(mat_double *in){
    int i,j;
    Rprintf("\n%dx%d matrix of doubles",in->n,in->p);
    for(i=0;i<in->n;i++){
	Rprintf("\n");
	for(j=0;j<in->p;j++)
	    /* printf("%d ", in->rows[i]->values[j]); */
	    Rprintf("%.3f ", mat_double_ij(in,i,j));
    }
    Rprintf("\n");
}






/* alternative print method for gsl vectors */
void print_gsl_vector(gsl_vector *in, char format[256]){
    int i;
    for(i=0;i<in->size;i++){
	Rprintf(format, in->data[i]);
    }
    Rprintf("\n");
}



/* check if an integer 'x' is in a vector of integers, and returns the matching position */
int in_vec_int(int x, vec_int *vec){
    int i=0;
    while(i<vec->length && x!=vec_int_i(vec, i)) i++; /* note: condition needs to be in this order */
    if(i==vec->length || vec->length<1) return -1; /* -1 will mean: no match*/
    return i;
}





/* permut the values of a vector of integers */
void permut_vec_int(vec_int *in, gsl_rng * rng){
    if(in->length<1) return;
    /* if(in->length != out->length){ */
    /* 	fprintf(stderr, "\n[in: matvec.c->permut_vec_int]\nInconsistent vector sizes: in = %d, out = %d",in->length,out->length); */
    /* 	exit(1); */
    /* } */

    gsl_ran_shuffle(rng, in->values, in->length, sizeof (int));
}





/*
   make one draw from a vector of probabilities
   returned value is (0,...,prob->length) with proba prob->values
*/
int draw_multinom(vec_double *prob, gsl_rng * rng){
    int i;
    double sumProb = sum_vec_double(prob), x, cumSum;

    /* draw 'x' between 0.0 and sumProb */
    x = gsl_rng_uniform(rng)*sumProb;

    /* get item id (id = i-1) */
    i=0;
    cumSum=0.0;
    do{
	cumSum += vec_double_i(prob,i++);
    } while(cumSum<x);

    return i-1;
} /* end draw_multinom */





/*
   same as above, using only the first 'n' items
*/
int draw_multinom_censored(vec_double *prob, int n, gsl_rng * rng){
    int i;
    double sumProb = 0.0, x, cumSum;

    for(i=0;i<n;i++){
	sumProb += vec_double_i(prob,i);
    }

    /* draw 'x' between 0.0 and sumProb */
    x = gsl_rng_uniform(rng)*sumProb;

    /* get item id (id = i-1) */
    i=0;
    cumSum=0.0;
    do{
	cumSum += vec_double_i(prob,i++);
    } while(cumSum<x);

    return i-1;
} /* end draw_multinom */






/* sample values of a vector of integers with/without replacement */
void sample_vec_int(vec_int *in, vec_int *out, bool replace, gsl_rng * rng){
    if(in->length<1 || out->length<1) return;
    if(out->length > in->length && !replace){
      error("\n[in: matvec.c->sample_vec_int]\nReplace is FALSE but sample size (%d) is bigger than input vector (%d)",out->length,in->length);
  	/* fprintf(stderr, "\n[in: matvec.c->sample_vec_int]\nReplace is FALSE but sample size (%d) is bigger than input vector (%d)",out->length,in->length); */
    	/* exit(1); */
    }

    if(replace){
	gsl_ran_sample(rng, out->values, out->length, in->values, in->length, sizeof (int));
    } else {
	gsl_ran_choose(rng, out->values, out->length, in->values, in->length, sizeof (int));
    }
} /* sample_vec_int */






/* draw values of a vector of integers with replacement using defined probabilities */
void draw_vec_int_multinom(vec_int *in, vec_int *out, vec_double *prob, gsl_rng * rng){
    int i;
    /* checks */
    if(in->length<1 || out->length<1) return;
    if(in->length != prob->length){
      error("\n[in: matvec.c->draw_vec_int_multinom]\nInput vector and vector of probabilities have different lengths (in:%d, prob:%d)", in->length, prob->length);
      /* fprintf(stderr, "\n[in: matvec.c->draw_vec_int_multinom]\nInput vector and vector of probabilities have different lengths (in:%d, prob:%d)", in->length, prob->length); */
      /* exit(1); */
    }

    /* draw values */
    for(i=0;i<out->length;i++){
	out->values[i] = vec_int_i(in, draw_multinom(prob, rng));
    }
    return;
} /* end draw_vec_int_multinom */





/* sort a vector of integers (ascending order)
   - in: input vector
   - out: vector of sorted values
   - idx: vector of indices (using R notation: out = in[idx])
*/
void sort_vec_int(vec_int *in, vec_int *out, vec_int *idx){
    if(out->length > in->length){
      error("\n[in: matvec.c->sort_vec_int]\nInput and output vectors have different lengths (in:%d out:%d)",in->length, out->length);
      /* fprintf(stderr, "\n[in: matvec.c->sort_vec_int]\nInput and output vectors have different lengths (in:%d out:%d)",in->length, out->length); */
      /* exit(1); */
    }

    int i, j, curMin, curMinIdx;
    idx->length=0;

    for(j=0;j<in->length;j++){
	/* printf("\n- sorting value %d\n",j); */
	/* find minimal value and its index, discarding already sorted indices */
	curMin=max_vec_int(in);
	curMinIdx=0;

	for(i=0;i<in->length;i++){
	    /* if(in_vec_int(i, idx)>=0){ */
	    /* 	printf("\nindex %d found in array idx (position:%d)", i, in_vec_int(i, idx)); */
	    /* } */
	    if(in_vec_int(i, idx)<0 && curMin>=vec_int_i(in,i)) {
		/* printf("\nentering the loop, i=%d\n",i); */
		/* printf("\nidx:"); print_vec_int(idx); */
		/* printf("\nmatch i in idx: %d\n", in_vec_int(i, idx)); */
		curMin=vec_int_i(in,i);
		curMinIdx = i;
	    }
	}

	/* printf("\n- minimum %d found at index %d",curMin,curMinIdx); */

	/* update vectors */
	out->values[j] = curMin;
	idx->values[j] = curMinIdx;
	idx->length = idx->length + 1;
    }

} /* end sort_vec_int */







/* sort a vector of doubles (ascending order)
   - in: input vector
   - out: vector of sorted values
   - idx: vector of indices (using R notation: out = in[idx])
*/
void sort_vec_double(vec_double *in, vec_double *out, vec_int *idx){
    if(out->length > in->length){
      error("\n[in: matvec.c->sort_vec_double]\nInput and output vectors have different lengths (in:%d out:%d)",in->length, out->length);
      /* fprintf(stderr, "\n[in: matvec.c->sort_vec_double]\nInput and output vectors have different lengths (in:%d out:%d)",in->length, out->length); */
      /* exit(1); */
    }

    int i, j, curMinIdx;
    double curMin;
    idx->length=0;

    for(j=0;j<in->length;j++){
	curMin=max_vec_double(in);
	curMinIdx=0;

	for(i=0;i<in->length;i++){
	    if(in_vec_int(i, idx)<0 && curMin>=vec_double_i(in,i)) {
		curMin=vec_double_i(in,i);
		curMinIdx = i;
	    }
	}

	/* update vectors */
	out->values[j] = curMin;
	idx->values[j] = curMinIdx;
	idx->length = idx->length + 1;
    }

} /* end sort_vec_double */









/*
   =========
   COPYING
   =========
*/

void copy_vec_int(vec_int *in, vec_int *out){
    int i;
    if(in->length != out->length){
      error("\n[in: matvec.c->copy_vec_int]\nInput and output vectors have different lengths (in:%d out:%d)",in->length, out->length);
      /* fprintf(stderr, "\n[in: matvec.c->copy_vec_int]\nInput and output vectors have different lengths (in:%d out:%d)",in->length, out->length); */
      /* exit(1); */
    }

    for(i=0;i<in->length;i++){
	out->values[i] = in->values[i];
    }
}



void copy_vec_double(vec_double *in, vec_double *out){
    int i;
    if(in->length != out->length){
      error("\n[in: matvec.c->copy_vec_double]\nInput and output vectors have different lengths (in:%d out:%d)",in->length, out->length);
      /* fprintf(stderr, "\n[in: matvec.c->copy_vec_double]\nInput and output vectors have different lengths (in:%d out:%d)",in->length, out->length); */
      /* exit(1); */
    }

    for(i=0;i<in->length;i++){
	out->values[i] = in->values[i];
    }
}




void copy_mat_int(mat_int *in, mat_int *out){
    int i;

    /* checks */
    if(in->n != out->n){
      error("\n[in: matvec.c->copy_mat_int]\nInput and output matrices have different numbers of rows (in:%d out:%d)",in->n, out->n);
      /* fprintf(stderr, "\n[in: matvec.c->copy_mat_int]\nInput and output matrices have different numbers of rows (in:%d out:%d)",in->n, out->n); */
      /* exit(1); */
    }
    if(in->p != out->p){
      error("\n[in: matvec.c->copy_mat_int]\nInput and output matrices have different numbers of columns (in:%d out:%d)",in->p, out->p);
      /* fprintf(stderr, "\n[in: matvec.c->copy_mat_int]\nInput and output matrices have different numbers of columns (in:%d out:%d)",in->p, out->p); */
      /* exit(1); */
    }

    /* copy */
    for(i=0;i<in->n;i++){
	copy_vec_int(in->rows[i],out->rows[i]);
    }
}




void copy_mat_double(mat_double *in, mat_double *out){
    int i;

    /* checks */
    if(in->n != out->n){
      error("\n[in: matvec.c->copy_mat_int]\nInput and output matrices have different numbers of rows (in:%d out:%d)",in->n, out->n);
      /* fprintf(stderr, "\n[in: matvec.c->copy_mat_int]\nInput and output matrices have different numbers of rows (in:%d out:%d)",in->n, out->n); */
      /* exit(1); */
    }
    if(in->p != out->p){
      error("\n[in: matvec.c->copy_mat_int]\nInput and output matrices have different numbers of columns (in:%d out:%d)",in->p, out->p);
	/* fprintf(stderr, "\n[in: matvec.c->copy_mat_int]\nInput and output matrices have different numbers of columns (in:%d out:%d)",in->p, out->p); */
    	/* exit(1); */
    }

    /* copy */
    for(i=0;i<in->n;i++){
	copy_vec_double(in->rows[i],out->rows[i]);
    }
}




/*
   ============
   BASIC STATS
   ============
*/

/*
   SUMS
*/
int sum_vec_int(vec_int *in){
    int i, out=0;
    for(i=0;i<in->length;i++){
	out += in->values[i];
    }
    return out;
}


double sum_vec_double(vec_double *in){
    int i;
    double out=0.0;
    for(i=0;i<in->length;i++){
	out += in->values[i];
    }
    return out;
}



/*
   MEANS
*/
double mean_vec_int(vec_int *in){
    double out=(double) sum_vec_int(in);
    return (out/in->length);
}


double mean_vec_double(vec_double *in){
    double out=(double) sum_vec_double(in);
    return (out/in->length);
}



/*
   MEDIANS
*/
double median_vec_int(vec_int *in){
    /* handle specific cases */
    if(in->length==0) return 0.0;
    if(in->length==1) return ((double) vec_int_i(in,0));

    /* allocate temporary vectors */
    double out=0.0;
    vec_int *sorted = alloc_vec_int(in->length), *idx=alloc_vec_int(in->length);

    /* sort values */
    sort_vec_int(in, sorted, idx);

    /* get median */
    out = 0.5 * (vec_int_i(in, floor(in->length/2.0)) + vec_int_i(in, 1+floor(in->length/2.0)));

    /* free and return */
    free_vec_int(idx);
    free_vec_int(sorted);
    return (out);
}


double median_vec_double(vec_double *in){
    /* handle specific cases */
    if(in->length==0) return 0.0;
    if(in->length==1) return ((double) vec_double_i(in,0));

    /* allocate temporary vectors */
    double out=0.0;
    vec_double *sorted = alloc_vec_double(in->length);
    vec_int *idx=alloc_vec_int(in->length);

    /* sort values */
    sort_vec_double(in, sorted, idx);

    /* get median */
    out = 0.5 * (vec_double_i(in, floor(in->length/2.0)) + vec_double_i(in, 1+floor(in->length/2.0)));

    /* free and return */
    free_vec_int(idx);
    free_vec_double(sorted);
    return (out);
}




/*
   MINIMUMS
*/

/* FIND MIN VALUE IN A VECTOR OF INTEGERS */
int min_vec_int(vec_int *vec){
    if(vec->length<1) return (int) NAN;
    int i, out=vec_int_i(vec,0);
    for(i=0;i<vec->length;i++) if(out>vec_int_i(vec,i)) out=vec_int_i(vec,i);
    return out;
}


/* FIND POSITION OF MIN VALUE IN A VECTOR OF INTEGERS */
int which_min_vec_int(vec_int *vec){
    if(vec->length<1) return (int) NAN;
    int i, curmin=vec_int_i(vec,0), out=0;
    for(i=0;i<vec->length;i++) {
	if(curmin>vec_int_i(vec,i)){
	    curmin=vec_int_i(vec,i);
	    out=i;
	}
    }
    return out;
}




/* FIND MIN VALUE IN A VECTOR OF DOUBLES */
double min_vec_double(vec_double *vec){
    if(vec->length<1) return (double) NAN;
    int i;
    double out=vec_double_i(vec,0);
    for(i=0;i<vec->length;i++) if(out>vec_double_i(vec,i)) out=vec_double_i(vec,i);
    return out;
}



/* FIND POSITION OF MIN VALUE IN A VECTOR OF DOUBLES */
int which_min_vec_double(vec_double *vec){
    if(vec->length<1) return (int) NAN;
    int i, out=0;
    double curmin=vec_double_i(vec,0);
    for(i=0;i<vec->length;i++) {
	if(curmin>vec_double_i(vec,i)){
	    curmin=vec_double_i(vec,i);
	    out=i;
	}
    }
    return out;
}






/*
   MAXIMUMS
*/

/* FIND MAX VALUE IN A VECTOR OF INTEGERS */
int max_vec_int(vec_int *vec){
    if(vec->length<1) return (int) NAN;
    int i, out=vec_int_i(vec,0);
    for(i=0;i<vec->length;i++) if(out<vec_int_i(vec,i)) out=vec_int_i(vec,i);
    return out;
}



/* FIND POSITION OF MAX VALUE IN A VECTOR OF INTEGERS */
int which_max_vec_int(vec_int *vec){
    if(vec->length<1) return (int) NAN;
    int i, curmax=vec_int_i(vec,0), out=0;
    for(i=0;i<vec->length;i++) {
	if(curmax<vec_int_i(vec,i)){
	    curmax=vec_int_i(vec,i);
	    out=i;
	}
    }
    return out;
}



/* FIND MAX VALUE IN A VECTOR OF DOUBLES */
double max_vec_double(vec_double *vec){
    if(vec->length<1) return (double) NAN;
    int i;
    double out=vec_double_i(vec,0);
    for(i=0;i<vec->length;i++) if(out<vec_double_i(vec,i)) out=vec_double_i(vec,i);
    return out;
}



/* FIND POSITION OF MAX VALUE IN A VECTOR OF DOUBLES */
int which_max_vec_double(vec_double *vec){
    if(vec->length<1) return (int) NAN;
    int i, out=0;
    double curmax=vec_double_i(vec,0);
    for(i=0;i<vec->length;i++) {
	if(curmax<vec_double_i(vec,i)){
	    curmax=vec_double_i(vec,i);
	    out=i;
	}
    }
    return out;
}







/* 
   ===================
   TRICKIER OPERATIONS
   ===================
*/

/* convolution for a positive discrete distribution
   (stored as a vector of doubles from 0 to t)
   a * b (t) = \sum_{u=0}^t a(t-u) b(u)
*/
void convol_vec_double(vec_double *in_a, vec_double *in_b, vec_double *out){
    int u, t;

    /* check sizes */
    if(in_a->length != in_b->length || in_a->length!=out->length){
      error("\n[in: matvec.c->convol_vec_double]\nInputs and output vectors have different lengths (in_a:%d in_b:%d out:%d)",in_a->length, in_b->length, out->length);
      /* fprintf(stderr, "\n[in: matvec.c->convol_vec_double]\nInputs and output vectors have different lengths (in_a:%d in_b:%d out:%d)",in_a->length, in_b->length, out->length); */
      /* exit(1); */
    }

    /* make computations */
    for(t=0;t<in_a->length;t++){
	for(u=0;u<=t;u++){
	    out->values[t] += vec_double_i(in_a, t-u) * vec_double_i(in_b, u);
	}
    }
} /* end convol_vec_double */





/*
  =========================
  === TESTING FUNCTIONS ===
  =========================
*/


/* int main(){ */
/*     /\* RANDOM NUMBER GENERATOR *\/ */
/*     time_t t = time(NULL); /\* time in seconds, used to change the seed of the random generator *\/ */
/*     const gsl_rng_type *typ; */
/*     gsl_rng_env_setup(); */
/*     typ=gsl_rng_default; */
/*     gsl_rng * rng=gsl_rng_alloc(typ); */
/*     gsl_rng_set(rng,t); /\* changes the seed of the random generator *\/ */

/*     int i, N = 5, P=10; */
/*     mat_int * test = alloc_mat_int(N, P); */

/*     print_mat_int(test); */
/*     free_mat_int(test); */

/*     vec_int *myVec = alloc_vec_int(30), *toto; */
/*     for(i=0;i<30;i++){ */
/* 	myVec->values[i] = 30-i; */
/*     } */
/*     printf("\nVector\n"); */
/*     print_vec_int(myVec); */
  
/*     vec_double *a = alloc_vec_double(4); */
/*     vec_double *convol = alloc_vec_double(4); */
/*     a->values[1]=3.0; */
/*     a->values[2]=2.0; */
/*     a->values[3]=4.0; */
/*     convol_vec_double(a,a,convol); */
/*     printf("\nvector 'a'\n"); */
/*     print_vec_double(a); */
/*     printf("\nConvolution (a*a)\n"); */
/*     print_vec_double(convol); */
 
/*     printf("\nMin/Max: %d, %d\n", min_vec_int(myVec), max_vec_int(myVec)); */
/*     printf("\nMin position/Max position: %d, %d\n", which_min_vec_int(myVec), which_max_vec_int(myVec)); */

/*     toto = alloc_vec_int(15); */
/*     sample_vec_int(myVec, toto, 1, rng); */
/*     printf("\n15 sampled values - with replacement \n"); */
/*     print_vec_int(toto); */

/*     sample_vec_int(myVec, toto, 1, rng); */
/*     printf("\nanother 15 sampled values - with replacement \n"); */
/*     print_vec_int(toto); */

/*     sample_vec_int(myVec, toto, 0, rng); */
/*     printf("\n15 sampled values - without replacement \n"); */
/*     print_vec_int(toto); */

/*     sample_vec_int(myVec, toto, 0, rng); */
/*     printf("\nanother 15 sampled values - without replacement \n"); */
/*     print_vec_int(toto); */

/*     permut_vec_int(myVec,rng); */
/*     printf("\npermut the vector myVec \n"); */
/*     print_vec_int(myVec); */

/*     permut_vec_int(myVec,rng); */
/*     printf("\nanother permutation of the vector myVec \n"); */
/*     print_vec_int(myVec); */
    
/*     printf("\n== sorting a vector ==\n"); */
/*     vec_int *idx, *sortedVec; */
/*     idx = alloc_vec_int(30); */
/*     vec_double *sorta = alloc_vec_double(a->length); */
/*     sortedVec = alloc_vec_int(30); */
/*     sort_vec_int(myVec, sortedVec, idx); */
/*     printf("\nvector to sort:"); */
/*     print_vec_int(myVec); */
/*     printf("\nsorted vector:"); */
/*     print_vec_int(sortedVec); */
/*     printf("\nindices:"); */
/*     print_vec_int(idx); */

/*     a->values[0]=300.012; */
/*     a->values[1]=-3.1; */
/*     a->values[2]=2.8; */
/*     a->values[3]=7.2; */
 
/*     printf("\nSorting a double\n");fflush(stdout); */
/*     printf("\nvector a:\n");fflush(stdout); */
/*     print_vec_double(a); */
/*     printf("\nsorted vector\n");fflush(stdout); */
/*     sort_vec_double(a, sorta, idx); */
/*     print_vec_double(sorta); */
/*     printf("\nmean of a: %f\n", mean_vec_double(a));fflush(stdout); */
/*     printf("\nmedian of a: %f\n", median_vec_double(a));fflush(stdout); */
    

/*     /\* copies *\/ */
/*     printf("\nCopy of sorted vector\n"); */
/*     vec_int *copyVec = alloc_vec_int(sortedVec->length); */
/*     copy_vec_int(sortedVec,copyVec); */
/*     print_vec_int(copyVec); */


/*     printf("\nCopy of a matrix\n"); */
/*     mat_double *mat = alloc_mat_double(2,2); */
/*     mat->rows[0]->values[0] = 1.1; */
/*     mat->rows[0]->values[1] = 2.1; */
/*     mat->rows[1]->values[1] = 666.0; */
/*     printf("\nmat\n"); */
/*     print_mat_double(mat); */
/*     mat_double *mat2 = alloc_mat_double(2,2); */
/*     printf("\nmat2 before copy\n"); */
/*     print_mat_double(mat2); */
/*     copy_mat_double(mat,mat2); */
/*     printf("\nmat2 after copy\n"); */
/*     print_mat_double(mat); */

/*     /\* vec_int *a = alloc_vec_int(10); *\/ */
/*     /\* for(i=0;i<10;i++){ *\/ */
/*     /\* 	a->values[i] = i; *\/ */
/*     /\* } *\/ */
/*     /\* printf("\nvector a: \n"); *\/ */
/*     /\* print_vec_int(a); *\/ */
/*     /\* for(i=0;i<10;i++){ *\/ */
/*     /\* 	printf("\n%d matches in a at position %d", i, in_vec_int(i,a)); *\/ */
/*     /\* } *\/ */

/*     /\* draw multinom *\/ */
/*     a->values[0] = 1.0; */
/*     a->values[1] = 1.0; */
/*     a->values[2] = 2.0; */
/*     a->values[3] = 4.0; */

/*     printf("\n\n30 draws from multinom prob=(1,1,2,4) \n");fflush(stdout); */
/*     for(i=0;i<30;i++){ */
/* 	printf("%d ", draw_multinom(a, rng));fflush(stdout); */
/*     } */

/*     printf("\n\n30 draws from multinom prob=(1,1,2,4) censored to n=3 items \n");fflush(stdout); */
/*     for(i=0;i<30;i++){ */
/* 	printf("%d ", draw_multinom_censored(a, 3, rng));fflush(stdout); */
/*     } */


/*     free_vec_int(toto); */
/*     free_vec_int(myVec); */
/*     free_vec_int(idx); */
/*     free_vec_int(sortedVec); */
/*     free_vec_int(copyVec); */
/*     free_vec_double(a); */
/*     free_vec_double(sorta); */
/*     free_vec_double(convol); */
/*     free_mat_double(mat); */
/*     free_mat_double(mat2); */
/*     gsl_rng_free(rng); */
/*     return 0; */
/* } */



/*
  gcc instructions

  gcc -o matvec matvec.c -lgsl -lgslcblas && ./matvec

  valgrind --leak-check=full matvec

*/

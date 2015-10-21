/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), January 2012.
  Licence: GPL >=2.

*/

#include "common.h"
#include "genclasses.h"

/* MACROS FOR DNAbin CLASS CONVERSION */
#define is_A(x)  (x==136)
#define is_G(x)  (x==72)
#define is_C(x)  (x==40)
#define is_T(x)  (x==24)





/*
  ============
  CONSTRUCTORS
  ============
*/

/* ALLOC DNASEQ OBJECT - ONE DNA SEQUENCE */
dnaseq * alloc_dnaseq(long int length){
    long int i;

    /* ALLOCATE OUTPUT */
    dnaseq *out = (dnaseq *) malloc(sizeof(dnaseq));
    if(out==NULL){
      error("\n[in: classes.c->alloc_dnaseq]\nNo memory left for creating DNA sequence. Exiting.\n");
      /* fprintf(stderr, "\n[in: classes.c->alloc_dnaseq]\nNo memory left for creating DNA sequence. Exiting.\n"); */
      /* exit(1); */
    }

    /* FILL/ALLOCATE CONTENT */
    out->seq =  (char*) malloc(length*sizeof(char));
    out->length = length;
    for(i=0;i<length;i++){
	out->seq[i] = '-';
    }

    return out;
}




/* ALLOC LIST_DNASEQ OBJECT - A LIST OF ALIGNED DNA SEQUENCES */
list_dnaseq * alloc_list_dnaseq(long int n, long int length){
    long int i;

    /* ALLOCATE OUTPUT */
    list_dnaseq *out = (list_dnaseq *) malloc(sizeof(list_dnaseq));
    if(out==NULL){
      error("\n[in: classes.c->alloc_list_dnaseq]\nNo memory left for creating list of DNA sequences. Exiting.\n");
      /* fprintf(stderr, "\n[in: classes.c->alloc_list_dnaseq]\nNo memory left for creating list of DNA sequences. Exiting.\n"); */
      /* exit(1); */
    }

    /* FILL/ALLOCATE CONTENT */
    out->list = (dnaseq **) malloc(n*sizeof(dnaseq *));
    if(out->list==NULL){
      error("\n[in: classes.c->alloc_list_dnaseq]\nNo memory left for creating list of DNA sequences. Exiting.\n");
      /* fprintf(stderr, "\n[in: classes.c->alloc_list_dnaseq]\nNo memory left for creating list of DNA sequences. Exiting.\n"); */
      /* exit(1); */
    }

    for(i=0;i<n;i++){
	out->list[i] = alloc_dnaseq(length);
    }
    out->n = n;
    out->length = length;
    return out;
}








/*
  ===========
  DESTRUCTORS
  ===========
*/

void free_dnaseq(dnaseq *in){
    if(in!=NULL){
	free(in->seq);
	free(in);
    }
}


void free_list_dnaseq(list_dnaseq *in){
    long int i;
    if(in!=NULL){
	for(i=0;i<in->n;i++){
	    free_dnaseq(in->list[i]);
	}
	free(in->list);
	free(in);
    }
}



/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/

void print_dnaseq(dnaseq *in){
    long int i;
    for(i=0;i<in->length;i++){
	Rprintf("%c", in->seq[i]);
    }
    Rprintf("\n");
}




void print_list_dnaseq(list_dnaseq *in){
    long int i;
    Rprintf("\nList of %d DNA sequences (size: %d)\n", in->n, in->length);
    for(i=0;i<in->n;i++){
	Rprintf("%d: ", i+1);
	print_dnaseq(in->list[i]);
    }
    Rprintf("\n");
}




/* convert a 'raw' of DNAbin object to char */
char DNAbin2char(unsigned char in){
    if(is_A(in)) return 'a';
    if(is_G(in)) return 'g';
    if(is_T(in)) return 't';
    if(is_C(in)) return 'c';
    return '-';
}



/*
  ==================
  EXTERNAL FUNCTIONS
  ==================
*/


/* convert DNAbin class to list_dnaseq */
list_dnaseq * DNAbin2list_dnaseq(unsigned char *in, long int *n, long int *length){
    long int i, j, count=0;

    /* ALLOC OUTPUT */
    list_dnaseq *out = alloc_list_dnaseq(*n,*length);

    /* FILL IN THE OUTPUT */
    for(i=0;i<*n;i++){
	for(j=0;j<*length;j++){
	    out->list[i]->seq[j] = DNAbin2char(in[count++]);
	}
    }

    /* printf("\nlist_dnaseq in C:\n"); */
    /* print_list_dnaseq(out); */
    /* printf("\n"); */

    /* RETURN */
    return out;
}




/* copy a dnaseq object */
void copy_dnaseq(dnaseq *in, dnaseq *out){
    if(out->length!=in->length){
      error("\n[in: genclasses.c->copy_dnaseq]\n.Input and output length mismatch (in:%d, out:%d)\n", in->length, out->length);
      /* fprintf(stderr, "\n[in: genclasses.c->copy_dnaseq]\n.Input and output length mismatch (in:%d, out:%d)\n", in->length, out->length); */
      /* exit(1); */
    }

    long int i;
    for(i=0;i<in->length;i++){
	out->seq[i] = in->seq[i];
    }
}




/*
  =======
  TESTING
  =======
*/


/* int main(){ */
/* 	const int N=5, L=30; */
/* 	int i,j; */

/* 	/\* alloc a list of sequences *\/ */
/* 	list_dnaseq * test = alloc_list_dnaseq(N, L); */

/* 	for(i=0;i<N;i++){ */
/* 		for(j=0;j<L;j++){ */
/* 			if(i*j % 5 ==0) test->list[i]->seq[j] = 'a'; */
/* 			else if(i*j % 3 ==0)test->list[i]->seq[j] = 't'; */
/* 			else if(i*j % 2 ==0)test->list[i]->seq[j] = 'g'; */
/* 			else test->list[i]->seq[j] = 'c'; */
/* 		} */
/* 	} */

/* 	print_list_dnaseq(test); */

/* 	/\* free and return *\/ */
/* 	free_list_dnaseq(test); */
/* 	return 0; */
/* } */


/* gcc instructions:

   gcc -o genclasses genclasses.c && ./genclasses


   valgrind --leak-check=full genclasses


*/




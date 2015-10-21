#include "common.h"
#include "matvec.h"
#include "genclasses.h"
#include "structures.h"


/*
  ======
   DATA
  ======
*/
data *alloc_data(int n, int nSeq, int length){
    /* allocate pointer */
    data *out = (data *) malloc(sizeof(data));
    if(out == NULL){
      error("\n[in: structures.c->alloc_data]\nNo memory left for creating data. Exiting.\n");
      /* fprintf(stderr, "\n[in: structures.c->alloc_data]\nNo memory left for creating data. Exiting.\n"); */
      /* exit(1); */
    }

    /* fill in integers */
    out->n = n;
    out->length = length;
    out->nSeq = nSeq;

    /* dates: collection times for each sequence */
    out->dates = alloc_vec_int(n);

    /* dna: list of DNA sequences */
    out->dna = alloc_list_dnaseq(nSeq, length);

    /* indices of DNA sequence for each case */
    out->idxCasesInDna = alloc_vec_int(n);

    /* locations of each case */
    out->locations = alloc_vec_int(n);

    return out;
} /* end alloc_data */





void free_data(data *in){
    free_vec_int(in->dates);
    free_list_dnaseq(in->dna);
    free_vec_int(in->idxCasesInDna);
    free_vec_int(in->locations);
    free(in);
} /* end free_data*/




void print_data(data *in){
    Rprintf("\n= Collection dates (timespan: %d)=\n",in->timespan);
    print_vec_int(in->dates);
    Rprintf("\n= Sequences =");
    print_list_dnaseq(in->dna);
    Rprintf("\n= Indices of DNA sequences for each case=\n");
    print_vec_int(in->idxCasesInDna);
    Rprintf("\n= Locations of the cases=\n");
    print_vec_int(in->locations);
} /* end print_data*/




/* Create a data object using inputs from R */
data * Rinput2data(unsigned char * DNAbinInput, int *Tcollec, int *n,
		   int *nSeq, int *length, int *idxCasesInDna, int *locations){
    int i, j, count=0;
    data * out = alloc_data(*n, *nSeq, *length);

    /* FILL IN VECTORS OF LENGTH N: DATES AND INDICES OF DNA */
    for(i=0;i<*n;i++){
	/* dates */
	out->dates->values[i] = Tcollec[i];

	/* indices of dna sequences for each case */
	/* -1 if no sequence (index is not in 0:(nSeq-1)) */
        out->idxCasesInDna->values[i] = (idxCasesInDna[i]<0||idxCasesInDna[i]>=*nSeq) ? -1 : idxCasesInDna[i];

	/* locations */
	out->locations->values[i] = locations[i];
    }

    out->timespan = max_vec_int(out->dates) - min_vec_int(out->dates);


    /* FILL IN DNA DATA */
    /* dna sequences */
    /* avoid using DNAbin2list_dnaseq here as it re-allocates memory */
    for(i=0;i<*nSeq;i++){
	for(j=0;j<*length;j++){
	    out->dna->list[i]->seq[j] = DNAbin2char(DNAbinInput[count++]);
	}
  }


    /* RETURN */
    return out;
} /* end Rinput2data */










/*
  =======
  GENTIME
  =======
*/

/* 'truncW' is the time at which w (generation time distrib.) is truncated to zero */
/* 'truncF' is the time at which f (time to collection distrib.) is truncated to zero */
/* 'maxK' is the maximum number of unobserved generations, for which we need convolutions */
/* 'maxT' must be larger than the largest time difference that can be observed between two related cases */
gentime *alloc_gentime(int maxK, int truncW, int truncF){
  /* allocate pointer */
    gentime *out = (gentime *) malloc(sizeof(gentime));
    if(out == NULL){
      error("\n[in: structures.c->alloc_gentime]\nNo memory left for creating gentime. Exiting.\n");
      /* fprintf(stderr, "\n[in: structures.c->alloc_gentime]\nNo memory left for creating gentime. Exiting.\n"); */
      /* exit(1); */
    }

    out->truncW = truncW>0 ? truncW : 1; /* make sur that p(0) is not zero */
    out->truncF = truncF>0 ? truncF : 1; /* make sur that p(0) is not zero */
    out->maxK = maxK>0 ? maxK : 1;

    /* allocate vector of densities */
    out->dens = alloc_mat_double(out->maxK, out->truncW*(out->maxK+2)); /* +2 to be on the safe side */
    out->collTime = alloc_vec_double(out->truncW*(out->maxK+2)); /* +2 to be on the safe side */

    /* return */
    return out;
} /* end alloc_gentime */




void free_gentime(gentime *in){
    free_mat_double(in->dens);
    free_vec_double(in->collTime);
    free(in);
} /* end free_gentime*/




void print_gentime(gentime *in){
    Rprintf("\n= Description of generation time function =\n");
    Rprintf("\n= Pre-computed density (truncated to 0 at %d)=\n",in->truncW);
    print_mat_double(in->dens);
    Rprintf("\n= Distribution of the time to collection =\n");
    print_vec_double(in->collTime);
} /* end print_gentime*/




/* get density from generation time funtion at time 't' with 'kappa_i' generations*/
double gentime_dens(gentime *in, int t, int kappa_i){
    /* error if requested kappa_i does not exist */
    if(kappa_i > in->maxK || kappa_i<1){
      error("\n[in: structures.c->gentime_dens]\nTrying to get density for %d generations (max: %d). Exiting.\n", kappa_i, in->maxK);
      /* fprintf(stderr, "\n[in: structures.c->gentime_dens]\nTrying to get density for %d generations (max: %d). Exiting.\n", kappa_i, in->maxK); */
      /* fflush(stdout); */
      /* exit(1); */
    }

    /* error if requested time too large */
    if(t >= in->maxK*in->truncW){
      error("\n[in: structures.c->gentime_dens]\nTrying to get density for %d time units (max: %d). Exiting.\n", t, in->maxK*in->truncW);
      /* fprintf(stderr, "\n[in: structures.c->gentime_dens]\nTrying to get density for %d time units (max: %d). Exiting.\n", t, in->maxK*in->truncW); */
      /* fflush(stdout); */
      /* exit(1); */
    }

    /* error if requested time is <= 0 */
    if(t <= 0){
      /* warning("\n[in: structures.c->gentime_dens]\nTrying to get density for %d time units (min: 1).\n", t); */
      return 0.0;
    }

    /* otherwise fetch density value */
    double out=mat_double_ij(in->dens, kappa_i-1, t);
    return out;
}




/* get density from the time to collection distribution */
double colltime_dens(gentime *in, int t){
    /* return 0 of t > truncF */
    if(t > in->truncF){
	return 0.0;
    }

    /* error if requested time is <= 0 */
    if(t <= 0){
      /* warning("\n[in: structures.c->colltime_dens]\nTrying to get density for %d time units (min: 1).\n", t); */
      return 0.0;
    }

    /* otherwise fetch density value */
    double out=vec_double_i(in->collTime, t);
    return out;
}



/*
 =======
  PARAM
 =======
*/

param *alloc_param(int n){
  /* allocate pointer */
    param *out = (param *) malloc(sizeof(param));
    if(out == NULL){
      error("\n[in: structures.c->alloc_param]\nNo memory left for creating param. Exiting.\n");
      /* fprintf(stderr, "\n[in: structures.c->alloc_param]\nNo memory left for creating param. Exiting.\n"); */
      /* exit(1); */
    }

    /* fill in integers */
    out->n = n;
    out->kappa_temp = 0;
    out->mut_model = 0;
    out->spa_model = 0;
    out->import_method = 0;

    /* allocates vectors of integers */
    out->Tinf = alloc_vec_int(n);
    out->alpha = alloc_vec_int(n);
    out->kappa = alloc_vec_int(n);

    /* fill in doubles */
    out->mu1 = 0.0001;
    out->mu1_prior = 0.0001;
    out->gamma = 1.0;
    out->pi = 1.0;
    out->pi_param1 = 0.0;
    out->pi_param2 = 0.0;
    out->spa_param1 = 0.0;
    out->spa_param2 = 0.0;
    out->spa_param1_prior = 0.0;
    out->spa_param2_prior = 0.0;
    out->outlier_threshold = 1000.0;
    out->phi = 0.5;
    out->phi_param1 = 1.0;
    out->phi_param2 = 1.0;

    /* return */
    return out;
} /* end alloc_param */




void free_param(param *in){
    free_vec_int(in->Tinf);
    free_vec_int(in->alpha);
    free_vec_int(in->kappa);
    free(in);
} /* end free_param*/




void print_param(param *in){
    Rprintf("\n= Tinf (infection dates) =\n");
    print_vec_int(in->Tinf);
    Rprintf("\n= Alpha_i (ancestries) =\n");
    print_vec_int(in->alpha);
    Rprintf("\n= Kappa_i (generations from nearest ancestor) =\n");
    print_vec_int(in->kappa);
    Rprintf("\n= mu1, mu2, gamma (mutation1, mutation2, coef, prior mu1) =\n");
    Rprintf("%.5f   %.5f   %.5f   %.5f", in->mu1, in->gamma*in->mu1, in->gamma, in->mu1_prior);
    Rprintf("\n= pi (proportion of observed cases) =\n");
    Rprintf("%.5f", in->pi);
    Rprintf("\n= priors on pi (parameter of beta distribution) =\n");
    Rprintf("%.5f  %.5f", in->pi_param1, in->pi_param2);
    Rprintf("\n= threshold used in imported case detection =\n");
    Rprintf("%.2f", in->outlier_threshold);
    Rprintf("\n= genetic model used =\n");
    Rprintf("%d", in->mut_model);
    Rprintf("\n= spatial model used (parameters) =\n");
    Rprintf("%d (param1=%.5f, param2=%.5f)", in->spa_model, in->spa_param1, in->spa_param2);
    Rprintf("\n= parameters of the spatial priors =\n");
    Rprintf("(mean prior param1=%.5f, mean prior param2=%.5f)", in->spa_param1_prior, in->spa_param2_prior);
    Rprintf("\n= imported case detection method used =\n");
    Rprintf("%d", in->import_method);
    Rprintf("\n= phi (proportion of nosocomial cases) =\n");
    Rprintf("%.5f", in->phi);
    Rprintf("\n= priors on phi (parameter of beta distribution) =\n");
    Rprintf("%.5f  %.5f", in->phi_param1, in->phi_param2);
} /* end print_param*/




void copy_param(param *in, param *out){
    /* copy atomic values */
    out->n = in->n;
    out->kappa_temp = in->kappa_temp;
    out->mu1 = in->mu1;
    out->mu1_prior = in->mu1_prior;
    out->gamma = in->gamma;
    out->pi = in->pi;
    out->pi_param1 = in->pi_param1;
    out->pi_param2 = in->pi_param2;
    out->phi = in->phi;
    out->phi_param1 = in->phi_param1;
    out->phi_param2 = in->phi_param2;
    out->spa_param1 = in->spa_param1;
    out->spa_param2 = in->spa_param2;
    out->spa_param1_prior = in->spa_param1_prior;
    out->spa_param2_prior = in->spa_param2_prior;
    out->outlier_threshold = in->outlier_threshold;
    out->mut_model = in->mut_model;
    out->spa_model = in->spa_model;
    out->import_method = in->import_method;

    copy_vec_int(in->Tinf,out->Tinf);
    copy_vec_int(in->alpha,out->alpha);
    copy_vec_int(in->kappa,out->kappa);
} /* end copy_param */









/*
 ============
  MCMC_PARAM
 ============
*/

mcmc_param *alloc_mcmc_param(int n){
  /* allocate pointer */
    mcmc_param *out = (mcmc_param *) malloc(sizeof(mcmc_param));
    if(out == NULL){
      error("\n[in: structures.c->alloc_mcmc_param]\nNo memory left for creating mcmc_param. Exiting.\n");
      /* fprintf(stderr, "\n[in: structures.c->alloc_mcmc_param]\nNo memory left for creating mcmc_param. Exiting.\n"); */
      /* exit(1); */
    }

    /* DETERMINE THE NUMBER OF Tinf */
    /* set to N/2, minimum 1 */
    out->n_move_Tinf = (int) n/2;
    out->n_move_Tinf = out->n_move_Tinf < 1 ? 1 : out->n_move_Tinf;


    /* DETERMINE THE NUMBER OF KAPPA AND ALPHA TO MOVE */
    /* set to N/2, minimum 1 */
    out->n_move_alpha = (int) n/2;
    out->n_move_alpha = out->n_move_alpha < 1 ? 1 : out->n_move_alpha;
    out->n_move_kappa = out->n_move_alpha;

    /* ALLOCATE VECTORS */
    out->idx_move_Tinf = alloc_vec_int(out->n_move_Tinf);
    out->idx_move_alpha = alloc_vec_int(out->n_move_alpha);
    out->idx_move_kappa = alloc_vec_int(out->n_move_kappa);
    out->all_idx = alloc_vec_int(n);
    out->candid_ances = alloc_vec_int(n+1);
    out->candid_ances_proba = alloc_vec_double(n+1);
    out->move_alpha = alloc_vec_double(n);
    out->move_kappa = alloc_vec_double(n);

    /* FILL IN INTEGERS */
    /* accept/reject counters */
    out->n_accept_mu1 = 0;
    out->n_reject_mu1 = 0;
    out->n_accept_gamma = 0;
    out->n_reject_gamma = 0;
    out->n_accept_spa1 = 0;
    out->n_reject_spa1 = 0;
    out->n_accept_spa2 = 0;
    out->n_reject_spa2 = 0;
    out->n_accept_Tinf = 0;
    out->n_reject_Tinf = 0;
    out->n_accept_alpha = 0;
    out->n_reject_alpha = 0;
    out->n_accept_kappa = 0;
    out->n_reject_kappa = 0;
    out->n_like_zero = 0;

    /* movement */
    out->move_mut = TRUE;
    out->move_pi = TRUE;
    out->move_phi = TRUE;
    out->move_spa = TRUE;
    /* out->move_phi = TRUE; */

    /* tuning */
    out->tune_any = TRUE;
    out->tune_mu1 = TRUE;
    out->tune_gamma = TRUE;
    out->tune_pi = TRUE;
    out->tune_phi = TRUE;
    out->step_notune = -1;

    /* misc */
    out->burnin=0;
    out->find_import_at=10000;
    out->find_import=TRUE;

    /* FILL IN DOUBLES */
    /* parameters of proposal */
    out->sigma_mu1 = 0.0;
    out->sigma_gamma = 0.0;
    out->sigma_pi = 0.0;
    out->sigma_phi = 0.0;
    out->sigma_spa1 = 0.0;
    out->sigma_spa2 = 0.0;
    out->lambda_Tinf = 0.0;


    /* RETURN */
    return out;
} /* end alloc_mcmc_param */




void free_mcmc_param(mcmc_param *in){
    free_vec_int(in->idx_move_Tinf);
    free_vec_int(in->idx_move_alpha);
    free_vec_int(in->idx_move_kappa);
    free_vec_int(in->all_idx);
    free_vec_int(in->candid_ances);
    free_vec_double(in->candid_ances_proba);
    free_vec_double(in->move_alpha);
    free_vec_double(in->move_kappa);
    free(in);
} /* end free_mcmc_param*/




void print_mcmc_param(mcmc_param *in){
    Rprintf("\nsigma for mu1: %.10f",in->sigma_mu1);
    Rprintf("\nsigma for gamma: %.10f",in->sigma_gamma);
    Rprintf("\nsigma for pi: %.10f",in->sigma_pi);
    Rprintf("\nsigma for phi: %.10f",in->sigma_phi);
    Rprintf("\nsigma for spa1: %.10f",in->sigma_spa1);
    Rprintf("\nsigma for spa2: %.10f",in->sigma_spa2);
    Rprintf("\nlambda for Tinf: %.10f",in->lambda_Tinf);
    Rprintf("\nnb moves for Tinf: %d",in->n_move_Tinf);
    Rprintf("\nnb moves for alpha: %d",in->n_move_alpha);
    Rprintf("\nnb moves for kappa: %d",in->n_move_kappa);

    Rprintf("\nmu1: nb. accepted: %d   nb. rejected: %d   (acc/rej ratio:%.3f)", in->n_accept_mu1, in->n_reject_mu1, (double) in->n_accept_mu1 / in->n_reject_mu1);

    Rprintf("\ngamma: nb. accepted: %d   nb. rejected: %d   (acc/rej ratio:%.3f)", in->n_accept_gamma, in->n_reject_gamma, (double) in->n_accept_gamma / in->n_reject_gamma);

    Rprintf("\npi: nb. accepted: %d   nb. rejected: %d   (acc/rej ratio:%.3f)", in->n_accept_pi, in->n_reject_pi, (double) in->n_accept_pi / in->n_reject_pi);

    Rprintf("\nphi: nb. accepted: %d   nb. rejected: %d   (acc/rej ratio:%.3f)", in->n_accept_phi, in->n_reject_phi, (double) in->n_accept_phi / in->n_reject_phi);

    Rprintf("\nspa1: nb. accepted: %d   nb. rejected: %d   (acc/rej ratio:%.3f)", in->n_accept_spa1, in->n_reject_spa1, (double) in->n_accept_spa1 / in->n_reject_spa1);

    Rprintf("\nspa2: nb. accepted: %d   nb. rejected: %d   (acc/rej ratio:%.3f)", in->n_accept_spa2, in->n_reject_spa2, (double) in->n_accept_spa2 / in->n_reject_spa2);

    Rprintf("\nTinf: nb. accepted: %d   nb. rejected: %d   (acc/rej ratio:%.3f)", in->n_accept_Tinf, in->n_reject_Tinf, (double) in->n_accept_Tinf / in->n_reject_Tinf);

    Rprintf("\nalpha: nb. accepted: %d   nb. rejected: %d   (acc/rej ratio:%.3f)", in->n_accept_alpha, in->n_reject_alpha, (double) in->n_accept_alpha / in->n_reject_alpha);

    Rprintf("\nkappa: nb. accepted: %d   nb. rejected: %d   (acc/rej ratio:%.3f)", in->n_accept_kappa, in->n_reject_kappa, (double) in->n_accept_kappa / in->n_reject_kappa);

    Rprintf("\nIndices of Tinf_i to move:\n");
    print_vec_int(in->idx_move_Tinf);
    Rprintf("\nIndices of alpha_i to move:\n");
    print_vec_int(in->idx_move_alpha);
    Rprintf("\nIndices of kappa_i to move:\n");
    print_vec_int(in->idx_move_kappa);

    Rprintf("\nVector of all indices (0:(n-1)):\n");
    print_vec_int(in->all_idx);

    Rprintf("\nVector of candidate ancestors:\n");
    print_vec_int(in->candid_ances);

    Rprintf("\nVector of candidate ancestors (proba):\n");
    print_vec_double(in->candid_ances_proba);

    Rprintf("\nTuned parameters: ");
    if(in->tune_mu1) Rprintf("mu1 ");
    if(in->tune_gamma) Rprintf("gamma ");
    if(in->tune_pi) Rprintf("pi ");
    if(in->tune_phi) Rprintf("phi ");
    if(in->tune_spa1) Rprintf("spa1 ");
    if(in->tune_spa2) Rprintf("spa2 ");
    Rprintf("\nTuning stopped at step %d\n", in->step_notune);

    Rprintf("\nMoved parameters:");
    if(in->move_mut) Rprintf("mu1 gamma ");
    /* if(in->move_alpha) Rprintf("alpha "); */
    /* if(in->move_kappa) Rprintf("kappa "); */
    if(in->move_Tinf) Rprintf("Tinf ");
    if(in->move_pi) Rprintf("pi ");
    if(in->move_phi) Rprintf("phi ");
    if(in->move_spa) Rprintf("spa ");
    Rprintf("\nMove alpha_i:");
    print_vec_double(in->move_alpha);
    Rprintf("\nMove kappa_i:");
    print_vec_double(in->move_kappa);
    if(in->find_import){
	Rprintf("\nFinding imported cases between chains %d and %d", in->burnin, in->find_import_at);
    }

} /* end print_mcmc_param */








/*
  Copy the content from one object to another.
  Allocation of the objects is done externally
*/
void copy_mcmc_param(mcmc_param *in, mcmc_param *out){
    /* copy atomic elements */
    out->n_accept_mu1 = in->n_accept_mu1;
    out->n_reject_mu1 = in->n_reject_mu1;
    out->n_accept_gamma = in->n_accept_gamma;
    out->n_reject_gamma = in->n_reject_gamma;
    out->n_accept_pi = in->n_accept_pi;
    out->n_reject_pi = in->n_reject_pi;
    out->n_accept_phi = in->n_accept_phi;
    out->n_reject_phi = in->n_reject_phi;
    out->n_accept_spa1 = in->n_accept_spa1;
    out->n_reject_spa1 = in->n_reject_spa1;
    out->n_accept_spa2 = in->n_accept_spa2;
    out->n_reject_spa2 = in->n_reject_spa2;
    out->n_accept_Tinf = in->n_accept_Tinf;
    out->n_reject_Tinf = in->n_reject_Tinf;
    out->n_accept_alpha = in->n_accept_alpha;
    out->n_reject_alpha = in->n_reject_alpha;
    out->n_accept_kappa = in->n_accept_kappa;
    out->n_reject_kappa = in->n_reject_kappa;

    out->n_move_Tinf = in->n_move_Tinf;
    out->n_move_alpha = in->n_move_alpha;
    out->n_move_kappa = in->n_move_kappa;

    out->sigma_mu1 = in->sigma_mu1;
    out->sigma_gamma = in->sigma_gamma;
    out->lambda_Tinf = in->lambda_Tinf;
    out->sigma_pi = in->sigma_pi;
    out->sigma_spa1 = in->sigma_spa1;
    out->sigma_spa2 = in->sigma_spa2;
    out->n_like_zero = in->n_like_zero;

    out->tune_any = in->tune_any;
    out->tune_mu1 = in->tune_mu1;
    out->tune_gamma = in->tune_gamma;
    out->tune_pi = in->tune_pi;
    out->tune_spa1 = in->tune_spa1;
    out->tune_spa2 = in->tune_spa2;
    out->step_notune = in->step_notune;

    out->move_mut = in->move_mut;
    out->move_pi = in->move_pi;
    out->move_phi = in->move_phi;
    out->move_spa = in->move_spa;
    out->burnin = in->burnin;
    out->find_import_at = in->find_import_at;
    out->find_import = in->find_import;

    /* copy vectors */
    copy_vec_int(in->idx_move_Tinf, out->idx_move_Tinf);
    copy_vec_int(in->idx_move_alpha, out->idx_move_alpha);
    copy_vec_int(in->idx_move_kappa, out->idx_move_kappa);
    copy_vec_int(in->all_idx, out->all_idx);
    copy_vec_int(in->candid_ances, out->candid_ances);
    copy_vec_double(in->candid_ances_proba, out->candid_ances_proba);
    copy_vec_double(in->move_alpha, out->move_alpha);
    copy_vec_double(in->move_kappa, out->move_kappa);
} /* end alloc_mcmc_param */






/*
  ======================
  >>>>> TESTING <<<<<
  ======================
*/



/* int main(){ */
/*     int i; */

/*     /\* data *\/ */
/*     data * dat = alloc_data(10,100); */
/*     printf("\nData\n"); */
/*     print_data(dat); */
/*     free_data(dat); */

/*     /\* gentime *\/ */
/*     gentime * gen = alloc_gentime(5, 20); */
/*     printf("\nGentime\n"); */
/*     print_gentime(gen); */
/*     printf("\nDensity for kappa=1, first 20 values:\n"); */
/*     for(i=0;i<20;i++) printf("%.6f ", gentime_dens(gen, i, 1)); */
/*     printf("\nDensity for kappa=2, first 20 values:\n"); */
/*     for(i=0;i<20;i++) printf("%.6f ", gentime_dens(gen, i, 2)); */
/*     printf("\nDensity for kappa=3, first 20 values:\n"); */
/*     for(i=0;i<20;i++) printf("%.6f ", gentime_dens(gen, i, 3)); */
/*     free_gentime(gen); */

/*     /\* param *\/ */
/*     param * par = alloc_param(10); */
/*     printf("\nParam\n"); */
/*     print_param(par); */
/*     free_param(par); */

/*     /\* /\\* mcmcParam *\\/ *\/ */
/*     mcmc_param * mcmcPar = alloc_mcmc_param(10); */
/*     printf("\nMcmcParam\n"); */
/*     print_mcmc_param(mcmcPar); */
/*     free_mcmc_param(mcmcPar); */


/*     printf("\n\n"); */
/*     return 0; */
/* } */



/*
   gcc instructions:

   gcc -o structures matvec.c genclasses.c structures.c -lgsl -lgslcblas -g -Wall

  ./structures

   valgrind --leak-check=full --track-origins=yes -v structures

*/

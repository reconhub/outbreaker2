#include "common.h"
#include "structures.h"
#include "genclasses.h"
#include "distances.h"
#include "init.h"
#include "prior.h"
#include "likelihood.h"
#include "moves.h"
#include "mcmc.h"




/*
  ======================
  MAIN EXTERNAL FUNCTION
  ======================
*/

void R_outbreaker(unsigned char *DNAbinInput, int *Tcollec, int *n, int *nSeq, int *length, 
		  int *idxCasesInDna, int *mutModel, double *gentimeDens, int *wTrunc, 
		  double *colltimeDens, int *fTrunc,
		  double *distMat, int *locations, int *spaModel,
		  int *ances, int *init_kappa, int *nIter, int *outputEvery, int *tuneEvery, 
		  double *piParam1, double *piParam2, 
		  double *phiParam1, double *phiParam2, 
		  double *initMu1, double *initGamma, 
		  double *initSpa1, double *initSpa2, 
		  double *spa1Prior, double *spa2Prior,
		  int *moveMut, int *moveAlpha, int *moveKappa, int *moveTinf, 
		  int *movePi, int *movePhi, int *moveSpa,
		  int *importMethod, int *findImportAt, int *burnin, 
		  double *outlierThreshold, int *maxK,
		  int *quiet, int *vecDist, int *stepStopTune,
		  char **resFileName, char **tuneFileName, int *seed){
    /* DECLARATIONS */
    int N = *n;
    gsl_rng *rng;
    data *dat;
    gentime *gen;
    param *par;
    dna_dist * dnaInfo;
    spatial_dist * spatialInfo;
    mcmc_param * mcmcPar;
    int i,j, counter;

    bool checkLike;
    bool findImport = (bool) *importMethod>0;

    /* INITIALIZE RNG */
    /* rng = create_gsl_rng((time_t) time(NULL)); */
    rng = create_gsl_rng((time_t) *seed);


    /* CONVERT DATA */
    dat = Rinput2data(DNAbinInput, Tcollec, n, nSeq, length, idxCasesInDna, locations);
    /* Rprintf("\n>>> Data <<<\n"); */
    /* print_data(dat); */


    /* /\* GET TIME SPAN *\/ */
    /* int TIMESPAN = max_vec_int(dat->dates) - min_vec_int(dat->dates); */
    /* Rprintf("\nTimespan is %d\n",TIMESPAN); */


    /* CREATE AND INIT GENERATION TIME */
    gen = alloc_gentime(*maxK, *wTrunc, *fTrunc);
    init_gentime(gen, gentimeDens, colltimeDens);
    /* Rprintf("\n>>> gentime info <<<\n"); */
    /* print_gentime(gen); */


    /* CREATE AND INIT PARAMETERS */
    par = alloc_param(N);
    init_param(par, dat,  gen, ances, init_kappa, *piParam1, *piParam2, *phiParam1, *phiParam2, *initMu1, *initGamma, *initSpa1, *initSpa2, *spa1Prior, *spa2Prior, *outlierThreshold, *mutModel, *spaModel, *importMethod, rng);
    /* Rprintf("\n>>> param <<<\n"); */
    /* print_param(par); */


    /* COMPUTE GENETIC DISTANCES */
    dnaInfo = compute_dna_distances(dat->dna, *mutModel);
    /* Rprintf("\n>>> DNA info <<<\n"); */
    /* print_dna_dist(dnaInfo); */


    /* CONVERT AND STORE SPATIAL DISTANCES */
    spatialInfo = doublevec2spatial_dist(distMat, n);
    /* Rprintf("\n>>> SPATIAL info <<<\n"); */
    /* print_spatial_dist(spatialInfo); */


    /* COMPUTE PRIORS */
    double logPrior = logprior_all(par);
    /* Rprintf("\nPrior value (log): %.10f\n", logPrior);/\* fflush(stdout); *\/ */

   /* COMPUTE LIKELIHOOD */
    double logLike = loglikelihood_all(dat, dnaInfo, spatialInfo, gen, par, rng);
    /* Rprintf("\n\n = Initial Log-likelihood value: %f\n", logLike); */

    /* COMPUTE POSTERIOR */
    double logPost = logposterior_all(dat, dnaInfo, spatialInfo, gen, par, rng);
    /* Rprintf("\nLog-posterior value: %.10f\n", logPost);/\* fflush(stdout); *\/ */

    /* ALLOCATE AND INITIALIZE MCMC PARAMETERS */
    /* Rprintf("\nBefore check init LL\n");/\* fflush(stdout); *\/ */

    mcmcPar = alloc_mcmc_param(N);
    init_mcmc_param(mcmcPar, par, dat, (bool) *moveMut, moveAlpha, moveKappa, (bool) *moveTinf, 
		    (bool) *movePi, (bool) *movePhi, (bool) *moveSpa, findImport, *burnin, *findImportAt);
    /* Rprintf("\nMCMC parameters\n");fflush(stdout); */
    /* print_mcmc_param(mcmcPar); */

    /* CHECK THAT INITIAL STATE HAS A NON-NULL LIKELIHOOD */
    checkLike = check_loglikelihood_all(dat, dnaInfo, spatialInfo, gen, par, rng);
    if(!checkLike){
      warning("\n\n!WARNING! Initial state of the chain has a likelihood of zero. The chain may never converge. Please consider using a different initial tree.\n");
    }

    /* Rprintf("\nAfter check init LL\n");/\* fflush(stdout); *\/ */
    /* Rprintf("\nBefore MCMC\n");/\* fflush(stdout); *\/ */

    /* RUN MCMC */
    mcmc(*nIter, *outputEvery, *resFileName, *tuneFileName, *tuneEvery,
	 (bool) *quiet, par, dat, dnaInfo, spatialInfo, gen, mcmcPar, rng);

    /* Rprintf("\nAfter MCMC\n");fflush(stdout); */

    /* FILL IN GENETIC DISTANCE VECTOR */
    counter = 0;
    for(i=0;i<(N-1);i++){
	for(j=i+1;j<N;j++){
	    vecDist[counter++] = mutation1_ij(i,j,dat,dnaInfo) + mutation2_ij(i,j,dat,dnaInfo);
	}
    }

    /* STORE STEP AT WHICH TUNING STOPPED */
    *stepStopTune = mcmcPar->step_notune;

    /* FREE MEMORY */
    free_data(dat);
    free_gentime(gen);
    free_param(par);
    free_dna_dist(dnaInfo);
    free_mcmc_param (mcmcPar);
    gsl_rng_free(rng);
} /* end R_outbreaker */





/* 

   Compilation instructions: 

   gcc -o outbreaker -Wall -g alloc.c matvec.c genclasses.c distances.c genlike.c logL.c prior.c moves.c mcmc.c init.c InputOutput.c tuneVariances.c outbreaker.c -lgsl -lgslcblas

   valgrind --leak-check=full outbreaker 

   valgrind -v --leak-check=full --track-origins=yes --show-reachable=yes outbreaker 


*/



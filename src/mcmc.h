
#ifndef __MCMC_H
#define __MCMC_H



/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/

void fprint_chains(FILE *file, data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, param *par, int step, gsl_rng *rng, bool quiet);

void fprint_mcmc_param(FILE *file, mcmc_param *mcmcPar, int step);




/*
   ================
   TUNING FUNCTIONS
   ================
*/

void tune_mu1(mcmc_param * in, gsl_rng *rng);

void tune_gamma(mcmc_param * in, gsl_rng *rng);

void tune_pi(mcmc_param * in, gsl_rng *rng);

void tune_phi(mcmc_param * in, gsl_rng *rng);

void tune_spa1(mcmc_param * in, gsl_rng *rng);

void tune_spa2(mcmc_param * in, gsl_rng *rng);


/* void tune_Tinf(mcmc_param * in, gsl_rng *rng); */




/*
   ===============================================
   METROPOLIS-HASTING ALGORITHM FOR ALL PARAMETERS
   ===============================================
*/
void mcmc_find_import(vec_int *areOutliers, int outEvery, int tuneEvery, bool quiet, param *par,
		      data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng);

void mcmc(int nIter, int outEvery, char outputFile[256], char mcmcOutputFile[256], int tuneEvery, bool quiet, param *par, data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng);


#endif

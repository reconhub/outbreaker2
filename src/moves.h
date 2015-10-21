
#ifndef __MOVES_H
#define __MOVES_H



/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/
bool look_for_aba(param *par, data *dat);

int choose_kappa_i(int T, gentime *gen, gsl_rng *rng);

int choose_alpha_i(int i, data *dat, param *currentPar, mcmc_param *mcmcPar, gsl_rng *rng);

int find_date_first_import(data *dat, param *par);

/*
  =====
  MOVES
  =====
*/

void move_mu1(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, mcmc_param *mcmcPar, gsl_rng *rng);

void move_gamma(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, mcmc_param *mcmcPar, gsl_rng *rng);

void move_pi(param *currentPar, param *tempPar, data *dat, mcmc_param *mcmcPar, gsl_rng *rng);

void move_phi(param *currentPar, param *tempPar, data *dat, spatial_dist *spaInfo, mcmc_param *mcmcPar, gsl_rng *rng);

void move_spa1(param *currentPar, param *tempPar, data *dat, spatial_dist *spaInfo, mcmc_param *mcmcPar, gsl_rng *rng);

void move_spa2(param *currentPar, param *tempPar, data *dat, spatial_dist *spaInfo, mcmc_param *mcmcPar, gsl_rng *rng);

void move_Tinf(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng);

void move_alpha_kappa(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng);

void move_Tinf_alpha_kappa(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng);

void swap_ancestries(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng);

/* NO LONGER USED */
/* void move_alpha(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng); */

/* void move_kappa(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng); */


#endif

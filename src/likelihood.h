#ifndef __LIKELIHOOD_H
#define __LIKELIHOOD_H



/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/
int find_sequenced_ancestor(int i, data *dat, dna_dist *dnaInfo, param *par);

int mutation1_ij(int i, int j, data *dat, dna_dist *dnaInfo);

int mutation2_ij(int i, int j, data *dat, dna_dist *dnaInfo);

int com_nucl_ij(int i, int j, data *dat, dna_dist *dnaInfo);

double gsl_ran_poisson_pdf_fixed(unsigned int k, double mu);

double proba_mut(int nbmut, int nbnucl, int kappa, double mu);

double log_proba_mut(int nbmut, int nbnucl, int kappa, double mu);


/*
  ====================
  LIKELIHOOD FUNCTIONS
  ====================
*/

double loglikelihood_i(int i, data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, param *par, gsl_rng *rng);

double loglikelihood_gen_i(int i, data *dat, dna_dist *dnaInfo, param *par, gsl_rng *rng);

double loglikelihood_spa_i(int i, data *dat, spatial_dist *spaInfo, param *par, gsl_rng *rng);

double loglikelihood_all(data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, param *par, gsl_rng *rng);

double loglikelihood_gen_all(data *dat, dna_dist *dnaInfo, param *par, gsl_rng *rng);

double loglikelihood_spa_all(data *dat, spatial_dist *spaInfo, param *par, gsl_rng *rng);

double loglike_kappa_all(param *par);

double loglikelihood_local_i(int i, data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, param *par, gsl_rng *rng);

/* double loglike_alpha_all(param *par); */

double logposterior_all(data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, param *par, gsl_rng *rng);

double sim_loglike_gen(data *dat, param *par, gsl_rng *rng);

bool check_loglikelihood_all(data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, param *par, gsl_rng *rng);

#endif

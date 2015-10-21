
/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk) and Anne Cori (a.cori@imperial.ac.uk), 2012.
  Licence: GPL >=2.
*/


/*
  This file contains the defition of the data structures used in this program, ad elementary methods (allocation, de-allocation, printing, accessing/setting content).
*/


#ifndef __STRUCTURES_H
#define __STRUCTURES_H

#include "common.h"
#include "matvec.h"
#include "genclasses.h"


/*
  =========
   CLASSES
  =========
*/

/* data: contains only observed data */
 typedef struct{
     int n, length, nSeq; /* n: number of observations; length: sequence length; nSeq: nb of sequences */
     vec_int * dates; /* collection dates*/
     list_dnaseq * dna; /* sequences */
     vec_int * idxCasesInDna; /* index of DNA sequence  in 'dna' for each case; -1 if DNA missing */
     int timespan; /* timespan of the data */
     vec_int * locations; /* integer indicating locations of cases - used in stratified dispersal model */
 } data;




/* descriptors/parameters of the generation time function 'w' and time to collection 'f' */
typedef struct{
    int truncW, truncF; /* value of truncation of the distributions*/
    int maxK; /* maximum value of kappa_i (i.e. max nb of generations between two cases) */
    mat_double *dens; /* pre-computed values of density: row 'i' gives densities for kappa=i at t=0, 1, ..., trunc-1*/
    vec_double *collTime; /* mass function of the "time interval to collection" */
} gentime;




/* param: contains augmented data, parameters and info about the model */
typedef struct{
    int n; /* number of observations, length of the vec_int objects */
    vec_int *Tinf; /* times of infection */
    vec_int *alpha; /* idx of closest ancestor for each case; -1 = unknown */
    vec_int *kappa; /* number of generations before a case and its closest ancestor */
    double mu1, mu1_prior; /* rate of transitions and its prior */
    double gamma; /* so that rate of transversions mu2 = gamma x mu1 */
    double pi; /* proportion of observed cases */
    double pi_param1, pi_param2; /* parameters of the Beta prior for pi */
    double spa_param1, spa_param2; /* parameters for the spatial model */
    double spa_param1_prior, spa_param2_prior; /* parameters of priors of spatial param */
    int kappa_temp; /* used to store temporary kappa for genetic LL */
    double outlier_threshold; /* threshold used in outlier detection */
    int mut_model; /* genetic model: 0=nothing; 1=1 mutation rate; 2=transi/transver */
    int spa_model; /* spatial model: 0=no, 1=exponential, 2=stratified exponential, ...*/
    int import_method; /* import method: 0=none; 1=based on genetic LL; 2=based on full LL */
    double phi; /* proba of nosocomial infection */
    double phi_param1, phi_param2; /* parameters of the Beta prior for phi */
} param;




/* mcmc_param: contains parameters of mcmc */
typedef struct{
    int n_accept_mu1, n_reject_mu1; /* accept/reject for mu1 */
    int n_accept_gamma, n_reject_gamma; /* accept/reject for gamma */
    int n_accept_pi, n_reject_pi; /* accept/reject for pi */
    int n_accept_phi, n_reject_phi; /* accept/reject for phi */
    int n_accept_spa1, n_reject_spa1; /* accept/reject for spa1 */
    int n_accept_spa2, n_reject_spa2; /* accept/reject for spa2 */
    int n_accept_Tinf, n_reject_Tinf; /* accept/reject for Tinf */
    int n_accept_alpha, n_reject_alpha; /* accept/reject for alpha */
    int n_accept_kappa, n_reject_kappa; /* accept/reject for kappa */
    int n_move_Tinf, n_move_alpha, n_move_kappa; /* number of Tinf, kappa and alpha to move at each chain */
    double sigma_mu1; /* sigma for normal proposal for mu1 */
    double sigma_gamma; /* sigma for normal proposal for gamma */
    double lambda_Tinf; /* lambda for Poisson movements of Tinf */
    double sigma_pi; /* sigma for normal proposal for pi */
    double sigma_phi; /* sigma for normal proposal for phi */
    double sigma_spa1, sigma_spa2; /* sigma for proposal for spa1 and spa2 */
    vec_int *idx_move_Tinf; /* vector of length n_move_Tinf giving indices of Tinf_i to move */
    vec_int *idx_move_alpha; /* vector of length n_move_alpha giving indices of alpha_i to move */
    vec_int *idx_move_kappa; /* vector of length n_move_kappa giving indices of kappa_i to move */
    vec_int *all_idx; /* vector of integers 0:(n-1) */
    vec_int *candid_ances; /* vector of candidate ancestors, used to move alpha_i */
    vec_double *candid_ances_proba; /* vector of proba for candidate ancestors, used to move alpha_i */
    /* double Pmove_alpha_old, Pmove_alpha_new; /\* used for accept ratio when moving alpha_i *\/ */
    int n_like_zero; /* number of times likelihood was zero */
    bool tune_any, tune_mu1, tune_gamma, tune_pi, tune_phi, tune_spa1, tune_spa2; /* logical indicating whether these proposals should be tuned */
    int step_notune; /* step at which all tuning stopped */
    bool move_mut, move_Tinf, move_pi, move_phi, move_spa; /* logical indicating what parameter should be moved */
    vec_double * move_alpha; /* vector indicating which alpha_i to move (1.0) or not (0.0) */
    vec_double * move_kappa; /* vector indicating which kappa_i to move (1.0) or not (0.0) */
    int burnin, find_import_at; /* chains between 'burnin' and 'find_import_at' are used to find imported cases */
    bool find_import; /* try to find and fix imported cases after 'burnin'? */
} mcmc_param;









/*
  =========
   METHODS
  =========
*/

/*
  ======
   DATA
  ======
*/
data *alloc_data(int n, int nSeq, int length);

void free_data(data *in);

void print_data(data *in);

data * Rinput2data(unsigned char * DNAbinInput, int *Tcollec, int *n,int *nSeq, int *length, int *idxCasesInDna, int *locations);




/*
  =======
  GENTIME
  =======
*/

gentime *alloc_gentime(int maxK, int truncW, int truncF);

void free_gentime(gentime *in);

void print_gentime(gentime *in);

double gentime_dens(gentime *in, int t, int kappa_i);

double colltime_dens(gentime *in, int t);


/*
 =======
  PARAM
 =======
*/

param *alloc_param(int n);

void free_param(param *in);

void print_param(param *in);

void copy_param(param *in, param *out);




/*
 ============
  MCMC_PARAM
 ============
*/

mcmc_param *alloc_mcmc_param(int n);

void free_mcmc_param(mcmc_param *in);

void print_mcmc_param(mcmc_param *in);

void copy_mcmc_param(mcmc_param *in, mcmc_param *out);








/* typedef struct{ */
/*     /\* to be estimated *\/ */
/*     gsl_matrix *beta; /\* person to person transmission rates between and within wards *\/ */
/*     double betaWardOut; /\* force of transmission from outside applied to patients in the wards *\/ */
/*     double betaOutOut; /\* force of transmission applied to individuals outside the wards *\/ */

/*     /\* double Sp; /\\* specificity of the test *\\/ assumed = 100% *\/ */
/*     double Se; /\* sensibility of the test *\/ */

/*     double Pi; /\* probability of being colonized at first admission *\/ */

/*     double mu; /\* mean duration of colonisation *\/ */
/*     double sigma; /\* std of the duration of colonisation *\/ */
/*     /\* *************** MAYBE NEED TO REPARAMETERIZE MU, SIGMA INTO MU, CV *************** *\/ */

/*     double nu1; /\* rate of transitions *\/ */
/*     /\* double nu2; /\\* rate of tranversions *\\/ *\/ */
/*     double kappa; /\* nu2 = kappa*nu1 *\/ */

/*     /\* double tau; /\\* time to the most recent common ancestor *\\/ *\/ */
/*     /\* double alpha; /\\* probability that two sampled isolates belong to the same lineage *\\/ *\/ */
/*     double weightNaGen; /\* weight used to replace missing genetic likelihood *\/ */
/* } parameters; */


/* typedef struct{ */
/*     double mu;  //mean duration of hospitalisation */
/*     double sigma;  //std of the duration of hospitalisation */
/* } hospDurationParam; */





/* typedef struct{ */
/*     long NbSimul; /\* nb of iterations of Metropolis-Hastings *\/ */
/*     int SubSample; /\* results recorded every SubSample iterations *\/ */
/*     int BurnIn; /\* nb of iterations considered as the Burn in period *\/ */

/*     /\* standard deviation for the proposition laws *\/ */
/*     gsl_matrix *Sigma_beta; */
/*     double Sigma_betaWardOut; */
/*     double Sigma_betaOutOut; */
/*     double Sigma_mu; */
/*     double Sigma_sigma; */
/*     double Sigma_nu1; */
/*     double Sigma_kappa; */
/*     double Sigma_tau; */
/*     double Sigma_alpha; */
/* } mcmcInternals; */





/* typedef struct{ */
/*     /\* average probabilities of acceptance *\/ */
/*     gsl_matrix *PourcAcc_beta; */
/*     double PourcAcc_betaWardOut; */
/*     double PourcAcc_betaOutOut; */
/*     double PourcAcc_mu; */
/*     double PourcAcc_sigma; */
/*     double PourcAcc_nu1; */
/*     double PourcAcc_kappa; */
/*     double PourcAcc_tau; */
/*     double PourcAcc_alpha; */

/* } acceptance; */





/* typedef struct{ */
/*     gsl_matrix *IsAccOK_beta; */
/*     double IsAccOK_betaWardOut; */
/*     double IsAccOK_betaOutOut; */
/*     double IsAccOK_mu; */
/*     double IsAccOK_sigma; */
/*     double IsAccOK_nu1; */
/*     double IsAccOK_kappa; */
/*     double IsAccOK_tau; */
/*     double IsAccOK_alpha; */
/* } isAcceptOK; */





/* typedef struct{ */
/*     gsl_matrix *NbProp_beta; */
/*     double NbProp_betaWardOut; */
/*     double NbProp_betaOutOut; */
/*     double NbProp_mu; */
/*     double NbProp_sigma; */
/*     double NbProp_nu1; */
/*     double NbProp_kappa; */
/*     double NbProp_tau; */
/*     double NbProp_alpha; */
/* } NbProposals; */





/* typedef struct{ */
/*     FILE *LogL; */
/*     FILE *Parameters; */
/*     FILE *ColonDates; */
/*     FILE *EndColonDates; */
/* } output_files; */



#endif

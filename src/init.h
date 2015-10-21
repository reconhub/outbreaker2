#ifndef __COMMON_H
#include "common.h"
#endif

#ifndef __INIT_H
#define __INIT_H

gsl_rng * create_gsl_rng(time_t t);

void init_gentime(gentime *in, double *Wvalues, double *Fvalues);

int find_maxLike_kappa_i(int T, gentime *gen);

void init_param(param *par, data *dat,  gentime *gen, int *ances, int *init_kappa, double pi_param1, double pi_param2, double phi_param1, double phi_param2, double init_mu1, double init_gamma, double init_spa1, double init_spa2, double spa1_prior, double spa2_prior, double outlier_threshold, int mut_model, int spa_model, int import_method, gsl_rng *rng);

void init_mcmc_param(mcmc_param *in, param *par, data *dat, bool move_mut, int *move_alpha, int *move_kappa, bool move_Tinf, bool move_pi, bool move_phi, bool move_spa, bool find_import, int burnin, int find_import_at);



#endif


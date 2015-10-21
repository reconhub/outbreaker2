
#ifndef __PRIOR_H
#define __PRIOR_H

void filter_logprob(double *in);

double gsl_ran_exponential_pdf_fixed(double x, double mu);

/* double logprior_alpha_i(int i, param *par); */

/* double logprior_kappa_i(int i, param *par); */

double logprior_pi(param *par);

double logprior_phi(param *par);

double logprior_spa1(param *par);

double logprior_spa2(param *par);

double logprior_mu1(param *par);

double logprior_gamma(param *par);

double logprior_all(param *par);


#endif

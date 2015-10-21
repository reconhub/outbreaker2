#include "common.h"
#include "structures.h"
#include "prior.h"


/* FUNCTION FILTERING LOG-PROB TO AVOID -INF or NaN*/
void filter_logprob(double *in){
    if(*in < NEARMINUSINF) *in = NEARMINUSINF;
}


/* fixed version of exponential density */
/* basic version returns NaN for mu=0; this one returns 0 */
double gsl_ran_exponential_pdf_fixed(double x, double mu){
    if(mu <= NEARZERO) return 0.0;
    return gsl_ran_exponential_pdf(x, mu);
}


/* /\* p(alpha_i) = ... */
/*    ...phi if alpha_i=0 */
/*    ...0 if alpha_i=i */
/*    ...(1-phi)/(n-1) in other cases */
/* *\/ */
/* double logprior_alpha_i(int i, param *par){ */
/*     /\* double out = (i==vec_int_i(par->alpha,i)) ? 0.0 : (double) 1.0/(par->n-1.0); *\/ */
/*     double out; */
/*     int ances = vec_int_i(par->alpha,i); */
/*     if(ances==0){ /\* external infection *\/ */
/* 	out = par->phi; */
/*     } else if(ances==i){ */
/* 	out = 0.0; /\* auto-infection (impossible) *\/ */
/*     } else { */
/* 	out =  (1.0 - par->phi)/(par->n - 1.0); /\* 'regular' infection *\/ */
/*     } */

/*     /\* put on log scale, filter, return *\/ */
/*     out = log(out); */
/*     filter_logprob(&out); */
/*     return out; */
/* } */



/* /\* p(kappa_i) = NB(1, 1-pi) with pi: prop obs cases *\/ */
/* double logprior_kappa_i(int i, param *par){ */
/*     double out = gsl_ran_negative_binomial_pdf((unsigned int) vec_int_i(par->kappa,i)-1, par->pi, 1.0); */
/*     out = log(out); */
/*     filter_logprob(&out); */
/*     return out; */
/* } */



/* p(mu_1) = exp(mu_init) */
/* two possible versions: unif or exponential */
/* exponential doesn't allow mu=0 */
double logprior_mu1(param *par){
    /* EXPONENTIAL VERSION */
    /* double out = gsl_ran_exponential_pdf_fixed(par->mu1, par->mu1_prior); */
    /* out = log(out); */
    /* filter_logprob(&out); */
    /* return out; */
    /* UNIFORM VERSION */
    return 0.0; /* log(1) = 0 */
}



/* p(gamma) = logN(1,1.25) */
double logprior_gamma(param *par){
    double out = gsl_ran_lognormal_pdf(par->gamma, 1.0, 1.25);
    out = log(out);
    filter_logprob(&out);
    return out;
}



/* p(pi) = beta(...,...) */
double logprior_pi(param *par){
    double out = gsl_ran_beta_pdf(par->pi, par->pi_param1, par->pi_param2);
    out = log(out);
    filter_logprob(&out);
    return out;
}



/* p(phi) = beta(...,...) */
double logprior_phi(param *par){
    double out = gsl_ran_beta_pdf(par->phi, par->phi_param1, par->phi_param2);
    out = log(out);
    filter_logprob(&out);
    return out;
}



/* prior on spatial param 1 */
double logprior_spa1(param *par){
    double out=0.0;

    /* SWITCH ACROSS MODELS */
    switch(par->spa_model){
	/* NULL MODEL - NO SPATIAL INFO */
    case 0:
	break;

	/* MODEL 1: exponential */
	/* (par->spa_param1 is the mean) */
	/* use an exponential prior for the mean of the exponential distribution */
    case 1:
	out = log(gsl_ran_exponential_pdf(par->spa_param1, par->spa_param1_prior));
	break;

	/* MODEL 2 */
    case 2:
      out = log(gsl_ran_exponential_pdf(par->spa_param1, par->spa_param1_prior));
      break;

	/* DEFAULT */
    default:
	break;
    }

    filter_logprob(&out);
    return out;
}



/* prior on spatial param 2 */
double logprior_spa2(param *par){
    double out=0.0;
    filter_logprob(&out);
    return out;
}



/* /\* p(phi) = beta(...,...) *\/ */
/* double logprior_phi(param *par){ */
/*     double out = gsl_ran_beta_pdf(par->phi, par->phi_param1, par->phi_param2); */
/*     out = log(out); */
/*     filter_logprob(&out); */
/*     return out; */
/* } */





double logprior_all(param *par){
    /* int i; */
    double out=0.0;

    out += logprior_mu1(par);
    if(par->mut_model>1) out += logprior_gamma(par);
    out += logprior_pi(par);
    if(par->spa_model>0){
	out += logprior_spa1(par);
	out += logprior_spa2(par);
    }
    if(par->spa_model==2){
	out += logprior_phi(par);
    }
    
    filter_logprob(&out);

    return(out);
}





/* int main(){ */
/*     double out; */
/*     /\* CREATE AND INIT PARAMETERS *\/ */
/*     param *par = alloc_param(3); */
    
/*     par->alpha->values[0] = -1; */
/*     par->alpha->values[1] = 0; */
/*     par->alpha->values[2] = 0; */
    
/*     par->kappa->values[0] = 1; */
/*     par->kappa->values[1] = 1; */
/*     par->kappa->values[2] = 1; */

/*     par->mu1 = 0.0001; */
/*     par->gamma = 1.0; */
    
/*     par->pi = 0.5; */

/*     printf("\nParameters\n"); */
/*     print_param(par); */

/*     /\* COMPUTE PRIORS *\/ */
/*     out = logprior_all(par); */
/*     printf("\nReturned value [log (prob)]: %.15f (%.15f)\n", out, exp(out)); */

/*     /\* FREE / RETURN *\/ */
/*     free_param(par); */
/*     return 0; */
/* } */



/*
  gcc instructions

  gcc -o prior matvec.c genclasses.c structures.c prior.c -lgsl -lgslcblas -Wall && ./prior

  valgrind --leak-check=full prior

*/


#include "common.h"
#include "matvec.h"
#include "genclasses.h"
#include "distances.h"
#include "structures.h"
#include "prior.h"
#include "likelihood.h"



/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/


/* FIND SEQUENCE TO USE FOR GENETIC LIKELIHOOD */
/* i: index of the case for which we seek a sequenced ancestor */
/* must return -1 if no ancestor sequenced (or i not sequenced) */
int find_sequenced_ancestor(int i, data *dat, dna_dist *dnaInfo, param *par){
  int nbNuclCommon = -1, curAnces = i, chainLength = 0;

    /* escape if no sequence for i */
    if(vec_int_i(dat->idxCasesInDna, i)<0) return -1;

    /* /\* debuging *\/ */
    /* printf("\nLooking for sequenced ancestor of %d, start with %d\n",i,vec_int_i(par->alpha,curAnces)); */
    /* fflush(stdout); */

    /* printf("%d",i);fflush(stdout); */

    /* store nb of generations from i to its closest sequenced ancestor */
    par->kappa_temp = 0;

    do{
	curAnces = vec_int_i(par->alpha,curAnces); /* move up the ancestry chain */
	/* printf("<-%d",curAnces);fflush(stdout); */
	if(curAnces>-1){
	  par->kappa_temp += vec_int_i(par->kappa,curAnces);
	  nbNuclCommon = com_nucl_ij(i, curAnces, dat, dnaInfo);
	}
	chainLength++;

	/* add warning if chain length exceeds n (likely loops) */
	if(chainLength > dat->n){
	  warning("\n\n!WARNING! Likely loops found in transmission chains when looking for sequenced ancestors.\n");
	}

	/* continue condition:
	   - sequenced ancestor found but no common nucleotide
	   - (and) current ancestor is not -1
	   - (and) we haven't gone too far back (indicative of a loop in transmission chain) */
    } while(nbNuclCommon<1 && curAnces>=0 && chainLength < dat->n);

   /* /\* debuging *\/ */
   /*  if(curAnces>0){ */
   /* 	printf("\nSequenced ancestor found for %d: %d (%d generations)\n",i,curAnces,par->kappa_temp); */
   /* 	fflush(stdout); */
   /*  } else { */
   /* 	printf("\nNo sequenced found for %d\n",i); */
   /* 	fflush(stdout); */
   /*  } */

    return curAnces;
} /* end find_sequenced_ancestor */





/* FIND NB TRANSITIONS BETWEEN CASES I AND J */
int mutation1_ij(int i, int j, data *dat, dna_dist *dnaInfo){
    /* if no nucleotide in common, return -1 */
    if(com_nucl_ij(i, j, dat, dnaInfo)<1) return -1;

    /* else read appropriate value in dnaInfo */
    return mat_int_ij(dnaInfo->mutation1, vec_int_i(dat->idxCasesInDna,i), vec_int_i(dat->idxCasesInDna,j));
} /* end mutation1_ij */





/* FIND NB TRANSVERSIONS BETWEEN CASES I AND J */
int mutation2_ij(int i, int j, data *dat, dna_dist *dnaInfo){
    /* if no nucleotide in common, return -1 */
    if(com_nucl_ij(i, j, dat, dnaInfo)<1) return -1;

    /* else read appropriate value in dnaInfo */
    return mat_int_ij(dnaInfo->mutation2, vec_int_i(dat->idxCasesInDna,i), vec_int_i(dat->idxCasesInDna,j));
} /* end mutation1_ij */






/* FIND NB OF COMPARABLE NUCLEOTIDES BETWEEN CASES I AND J */
int com_nucl_ij(int i, int j, data *dat, dna_dist *dnaInfo){
    /* return -1 if i or j is unknown case */
    if(i<0 || j<0) return -1;

    /* if 1 missing sequence, return -1 */
    if(vec_int_i(dat->idxCasesInDna,i)<0 || vec_int_i(dat->idxCasesInDna,j)<0) return -1;

    return mat_int_ij(dnaInfo->nbcommon, vec_int_i(dat->idxCasesInDna,i), vec_int_i(dat->idxCasesInDna,j));
} /* end mutation1_ij */



/* Fixed version of gsl_ran_poisson_pdf */
/* Original returns P(0|0) = NaN; this one returns 1.0 */
double gsl_ran_poisson_pdf_fixed(unsigned int k, double mu){
    if(mu <= NEARZERO && k==0) return 1.0;
    return gsl_ran_poisson_pdf(k, mu);
} /* gsl_ran_poisson_pdf_fixed */




/* Compute proba of ... mutations */
/* - binomial density without the binomial coefficient - */
/* Assume that no site can mutate twice over kappa generations */
/* there are kappa*nbnucl 'draws' (possible mutation occurences) */
/* Formula: p = mu^nbmut x (1-mu)^{(kappa x nbnucl) - nbmut} */
double proba_mut(int nbmut, int nbnucl, int kappa, double mu){
    double out=0.0;
    /* old version, numerical approx problems for nbmut high */
    out = gsl_sf_pow_int(mu, nbmut) * gsl_sf_pow_int(1.0-mu, kappa*nbnucl - nbmut);
    return out;
}

/* Compute log-proba of ... mutations */
/* - binomial density without the binomial coefficient - */
/* Assume that no site can mutate twice over kappa generations */
/* there are kappa*nbnucl 'draws' (possible mutation occurences) */
/* Formula: p = mu^nbmut x (1-mu)^{(kappa x nbnucl) - nbmut} */
double log_proba_mut(int nbmut, int nbnucl, int kappa, double mu){
    double out=0.0;
    /* note: integers are automatically promoted to double here */
    out = nbmut * log(mu) + (kappa*nbnucl - nbmut) * log(1.0-mu);
    return out;
}





/*
  ====================
  LIKELIHOOD FUNCTIONS
  ====================
*/

/* LOG-LIKELIHOOD FOR INDIVIDUAL 'i' */
double loglikelihood_i(int i, data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, param *par, gsl_rng *rng){
    int ances=vec_int_i(par->alpha,i);
    double out=0.0;


    /* = EXTERNAL CASES = */
    if(ances < 0){
      /* PROBA OF SAMPLING TIME */
      if(vec_int_i(dat->dates,i) <= vec_int_i(par->Tinf,i)){ /* fool proof */
	out += NEARMINUSINF;
      } else {
	out = log(colltime_dens(gen, vec_int_i(dat->dates,i) - vec_int_i(par->Tinf,i)));
      }

      /* /\* PROBA OF EXTERNAL CASE *\/ */
      /* out += log(par->phi); */

      /* PROBA OF INFECTION TIME (UNIFORM OVER TIMESPAN) */
      out -= log((double) dat->timespan);

      /* FILTER AND RETURN */
      filter_logprob(&out);
      /* printf("\nlikelihood (imported case): %f\n", out);fflush(stdout); */
      return out;
    }


    /* = INTERNAL CASES = */
    /* GENETIC LIKELIHOOD */
    out += loglikelihood_gen_i(i, dat, dnaInfo, par, rng);

    /* EPIDEMIOLOGICAL LIKELIHOOD */
    /* LIKELIHOOD OF COLLECTION DATE */
    if(vec_int_i(dat->dates,i) <= vec_int_i(par->Tinf,i)){ /* fool proof */
      out += NEARMINUSINF;
    } else {
      out += log(colltime_dens(gen, vec_int_i(dat->dates,i) - vec_int_i(par->Tinf,i)));
    }

    /* LIKELIHOOD OF INFECTION TIME */
    /* printf("\ninfection date: %.10f\n", log(gentime_dens(gen, vec_int_i(par->Tinf,i) - vec_int_i(par->Tinf,ances), vec_int_i(par->kappa,i)))); */
    if(vec_int_i(par->Tinf,i) <= vec_int_i(par->Tinf,ances)){ /* fool proof */
      out += NEARMINUSINF;
    } else {
      out += log(gentime_dens(gen, vec_int_i(par->Tinf,i) - vec_int_i(par->Tinf,ances), vec_int_i(par->kappa,i)));
    }

    /* PROBA OF (KAPPA_I-1) UNOBSERVED CASES */
    if(vec_int_i(par->kappa,i)<1 || par->pi<=0 || par->pi >1){ /* fool proof */
      out += NEARMINUSINF;
    } else {
      out += log(gsl_ran_negative_binomial_pdf((unsigned int) vec_int_i(par->kappa,i)-1, par->pi, 1.0));
    }

    /* SPATIAL LIKELIHOOD */
    out += loglikelihood_spa_i(i, dat, spaInfo, par, rng);

    /* FILTER AND RETURN */
    filter_logprob(&out);

    /* printf("\nlikelihood (internal case): %f\n", out);fflush(stdout); */

    return out;
} /* end loglikelihood_i */





/* GENETIC LOG-LIKELIHOOD FOR INDIVIDUAL 'i' */
double loglikelihood_gen_i(int i, data *dat, dna_dist *dnaInfo, param *par, gsl_rng *rng){
    int ances;
    double out=0.0;

    /* ESCAPE OF NO EVOLUTION MODEL CHOSEN */
    if(par->mut_model==0) return 0.0;

    /* ESCAPE IF NOT SEQUENCE FOR I */
    if(vec_int_i(dat->idxCasesInDna, i)<0) return 0.0;

    /* FIND MOST RECENT SEQUENCED ANCESTOR */
    /* Rprintf("\nSeeking most recent sequenced ancestor of %d",i); */
    ances = find_sequenced_ancestor(i, dat, dnaInfo, par);
    /* Rprintf("\n...found %d",ances); */

    /* NO DNA INFO AVAIL (IMPORTED CASES/MISSING SEQUENCES) */
    if(ances < 0) {
	/* return sim_loglike_gen(dat, par, rng); */
	return 0.0;
    }

    /* SWITCH ACROSS MODELS */
    switch(par->mut_model){
	/* MODEL 1: only one type of mutations */
    case 1:
      if(com_nucl_ij(i, ances, dat, dnaInfo)>0){
	if(mutation1_ij(i, ances, dat, dnaInfo)<0 || par->kappa_temp<0 || par->mu1<0 || par->mu1 >1){ /* fool proof */
	  out += NEARMINUSINF;
	} else {
	  out += log_proba_mut(mutation1_ij(i, ances, dat, dnaInfo), com_nucl_ij(i, ances, dat, dnaInfo), par->kappa_temp, par->mu1);
	}
      }
      break;

  /* MODEL 2: transitions and transversions */
    case 2:
      if(com_nucl_ij(i, ances, dat, dnaInfo)>0){
	if(mutation1_ij(i, ances, dat, dnaInfo)<0 || mutation2_ij(i, ances, dat, dnaInfo)<0 || par->kappa_temp<0 || par->mu1<0 || par->mu1 >1 || par->gamma<0){ /* fool proof */
	  out += NEARMINUSINF;
	} else {
	  /* transitions */
	  out += log_proba_mut(mutation1_ij(i, ances, dat, dnaInfo), com_nucl_ij(i, ances, dat, dnaInfo), par->kappa_temp, par->mu1);

	  /* transversions */
	  out += log_proba_mut(mutation2_ij(i, ances, dat, dnaInfo), com_nucl_ij(i, ances, dat, dnaInfo), par->kappa_temp,  par->gamma * par->mu1);
	}
      }
      break;
	/* DEFAULT */
    default:
	break;
    } /* end switch */

    /* RETURN */
    filter_logprob(&out);

    return out;
} /* end loglikelihood_gen_i */





/* SPATIAL LOG-LIKELIHOOD FOR INDIVIDUAL 'i' */
double loglikelihood_spa_i(int i, data *dat, spatial_dist *spaInfo, param *par, gsl_rng *rng){
    double out=0.0, Dij=0.0;
    int ances=0;

    /* SWITCH ACROSS MODELS */
    switch(par->spa_model){
	/* NULL MODEL - NO SPATIAL INFO */
    case 0:
	break;

	/* MODEL 1: exponential */
	/* f(x, rate) = rate e^{-rate*x} */
	/* mean = 1/rate */
	/* par->spa_param1 is the mean */
    case 1:
	/* printf("\nLooking for spa like for %d\n", i);fflush(stdout); */
	ances = vec_int_i(par->alpha, i);
	if(ances>=0){ /* only if not imported case */
	  /* printf("\nancestor: %d\n", ances);fflush(stdout); */
	  Dij = get_spatial_dist(spaInfo, ances, i);
	  if(Dij < 0 || par->spa_param1<0) { /* fool-proof */
	    out = NEARMINUSINF;
	  } else {
	    /* printf("\nDistance: %.5f\n", Dij);fflush(stdout); */
	    out = log(gsl_ran_exponential_pdf(Dij, par->spa_param1));
	    /* printf("\nLog-like: %.5f\n", out);fflush(stdout); */
	  }
	}
	break;

	/* MODEL 2: stratified exponential */
	/* Either local/nosocomial transmission, or exponential diffusion */
	/* phi: proba of local/nosocomial transmisson */
    /* case 2: */
    /* 	ances = vec_int_i(par->alpha, i); */
    /* 	if(ances>=0){ /\* only if not imported case *\/ */
    /* 	    /\* if same location / nosocomial transmission *\/ */
    /* 	    if(vec_int_i(dat->locations,i)==vec_int_i(dat->locations,ances)){ */
    /* 	      /\* printf("\nlocations of %d and %d are the same", i, ances);fflush(stdout); *\/ */
    /* 	      /\* printf("\nvalue of phi:%.5f", par->phi);fflush(stdout); *\/ */
    /* 	      out = log(par->phi); */
    /* 	    } else { */
    /* 		/\* different locations - non nosocomial *\/ */
    /* 		Dij = get_spatial_dist(spaInfo, ances, i); */
    /* 		out = log((1.0 - par->phi) * gsl_ran_exponential_pdf(Dij, par->spa_param1)); */
    /* 		/\* printf("\nlocations of %d and %d are different", i, ances);fflush(stdout); *\/ */
    /* 		/\* printf("\nvalue of LL:%.5f", out);fflush(stdout); *\/ */
    /* 		/\* printf("\ndetail:(1.0-phi)=%.5f ; p(Dij)=%.15f ; spa1=%.5f", 1.0 - par->phi, gsl_ran_exponential_pdf(Dij, par->spa_param1), par->spa_param1);fflush(stdout); *\/ */
    /* 	    } */
    /* 	} */
    /* 	break; */

	/* DEFAULT */
    default:
	break;
    }

    /* RETURN */
    filter_logprob(&out);

    return out;
} /* end loglikelihood_spa_i*/





/* LOG-LIKELIHOOD FOR ALL INDIVIDUALS */
double loglikelihood_all(data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, param *par, gsl_rng *rng){
    int i;
    double out=0.0;

    for(i=0;i<dat->n;i++){
	out += loglikelihood_i(i, dat, dnaInfo, spaInfo, gen, par, rng);
    }

    filter_logprob(&out);

    return out;
}







/* GENETIC LOG-LIKELIHOOD FOR ALL INDIVIDUALS */
double loglikelihood_gen_all(data *dat, dna_dist *dnaInfo, param *par, gsl_rng *rng){
    int i;
    double out=0.0;

    /* ESCAPE OF NO EVOLUTION MODEL CHOSEN */
    if(par->mut_model==0) return 0.0;

    for(i=0;i<dat->n;i++){
	out += loglikelihood_gen_i(i, dat, dnaInfo, par, rng);
    }

    filter_logprob(&out);

    return out;
}






/* SPATIAL LOG-LIKELIHOOD FOR ALL INDIVIDUALS */
double loglikelihood_spa_all(data *dat, spatial_dist *spaInfo, param *par, gsl_rng *rng){
    int i;
    double out=0.0;

    for(i=0;i<dat->n;i++){
	out += loglikelihood_spa_i(i, dat, spaInfo, par, rng);
    }

    filter_logprob(&out);

    return out;
}






/* LOG-LIKELIHOOD KAPPA FOR ALL INDIVIDUALS */
double loglike_kappa_all(param *par){
    double out = 0.0;
    int i;
    for(i=0;i<par->n;i++){
      if(vec_int_i(par->kappa,i)<1){ /* fool-proof */
	out += NEARMINUSINF;
      } else {
	out += log(gsl_ran_negative_binomial_pdf((unsigned int) vec_int_i(par->kappa,i)-1, par->pi, 1.0));
      }
    }
    filter_logprob(&out);
    return out;
}






/* LOCAL LOG-LIKELIHOOD FOR INDIVIDUAL I */
/* only computed on elements changed by a movement of i */
/* - p(i|ancestor) */
/* - p(i's descendents|i) */
double loglikelihood_local_i(int i, data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, param *par, gsl_rng *rng){
    int j;
    double out=0.0;

    /* i's likelihood */
    out += loglikelihood_i(i, dat, dnaInfo, spaInfo, gen, par, rng);

    for(j=0;j<dat->n;j++){
	/* likelihoods of i's descendents */
	if(i!=j && vec_int_i(par->alpha, j)==i) out += loglikelihood_i(j, dat, dnaInfo, spaInfo, gen, par, rng);
    }

    filter_logprob(&out);

    return out;
}


/* /\* LOG-LIKELIHOOD ALPHA FOR ALL INDIVIDUALS *\/ */
/* double loglike_alpha_all(param *par){ */
/*     double out = 0.0; */
/*     int i; */
/*     for(i=0;i<par->n;i++){ */
/* 	out += (vec_int_i(par->alpha,i)<0) ? log(par->phi) : log(1.0-par->phi); */
/*     } */
/*     filter_logprob(&out); */
/*     return out; */
/* } */





/* LOG-POSTERIOR FOR ALL INDIVIDUALS */
double logposterior_all(data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, param *par, gsl_rng *rng){
    double out = logprior_all(par) + loglikelihood_all(dat, dnaInfo, spaInfo, gen, par, rng);

    filter_logprob(&out);

    return out;
}





/* SIMULATE GENETIC LOG-LIKELIHOOD */
double sim_loglike_gen(data *dat, param *par, gsl_rng *rng){
    double out = 0.0, lambda1, lambda2;

    /* compute param of the Poisson distributions */
    lambda1 = (double) dat->length * par->mu1;
    lambda2 = (double) dat->length * par->mu1 * par->gamma;

    /* compute log-likelihood */
    out += log(gsl_ran_poisson_pdf_fixed(gsl_ran_poisson(rng, lambda1) , lambda1));
    out += log(gsl_ran_poisson_pdf_fixed(gsl_ran_poisson(rng, lambda2) , lambda2));

    /* /\* penalize the likelihoood *\/ */
    /* out -= 10; */

    /* filter and return */
    filter_logprob(&out);

    return out;
} /* end sim_loglike_gen */





/* CHECK LOG-LIKELIHOOD FOR ALL INDIVIDUALS */
/* returns TRUE if all is fine, FALSE if likelihood is zero */
bool check_loglikelihood_all(data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, param *par, gsl_rng *rng){
    int i, ances;
    double temp;
    bool out=TRUE;

    for(i=0;i<dat->n;i++){
	temp = loglikelihood_i(i, dat, dnaInfo, spaInfo, gen, par, rng);
	filter_logprob(&temp);

	if(temp <= NEARMINUSINF){
	    out = FALSE;
	    Rprintf("\nlikelihood for ancestry of %d is zero", i+1);
	    /* fflush(stdout); */

	    /* display genetic likelihood */
	    temp = loglikelihood_gen_i(i, dat, dnaInfo, par, rng);
	    filter_logprob(&temp);
	    Rprintf("\ni=%d: genetic log-like is: %f", i+1, temp);
	    /* fflush(stdout); */
	    if(temp <= NEARMINUSINF) Rprintf(" (i.e., zero)");
	    /* fflush(stdout); */

	    /* display epi likelihood */
	    ances=vec_int_i(par->alpha,i);

	    /* likelihood of collection date */
	    temp = log(colltime_dens(gen, vec_int_i(dat->dates,i) - vec_int_i(par->Tinf,i)));
	    filter_logprob(&temp);
	    Rprintf("\ni=%d: collection date (t_%d=%d,Tinf_%d=%d) log-like is: %f", i+1, i+1, vec_int_i(dat->dates,i), i+1, vec_int_i(par->Tinf,i), temp);
	    /* fflush(stdout); */
	    if(temp <= NEARMINUSINF) Rprintf(" (i.e., zero)");
	    /* fflush(stdout); */

	    /* likelihood of infection time */
	    temp = log(gentime_dens(gen, vec_int_i(par->Tinf,i) - vec_int_i(par->Tinf,ances), vec_int_i(par->kappa,i)));
	    filter_logprob(&temp);
	    Rprintf("\ni=%d: infection time (Tinf=%d,Tances=%d) log-like is: %f", i+1, vec_int_i(par->Tinf,i),vec_int_i(par->Tinf,ances), temp);
	    /* fflush(stdout); */

	    /* spatial likelihood*/
	    temp = loglikelihood_spa_i(i, dat, spaInfo, par, rng);
	    filter_logprob(&temp);
	    Rprintf("\ni=%d: spatial log-like is: %f", i+1, temp);

	    if(temp <= NEARMINUSINF) Rprintf(" (i.e., zero)");
	    /* fflush(stdout); */

	}
    }

    return out;
} /* end check_loglikelihood_all */








/*
>>>> TESTING <<<<
*/

/* #include "init.h" */
/* int main(){ */
/*   /\* DECLARATIONS *\/ */
/*     int TIMESPAN; */
/*     data *dat; */
/*     gentime *gen; */
/*     param *par; */
/*     dna_dist * dnaInfo; */

/*     double logPrior, logLike, logPost; */



/*     /\* CONVERT DATA *\/ */
/*     dat = alloc_data(3,10); */
/*     dat->dates->values[0] = 0; */
/*     dat->dates->values[1] = 2; */
/*     dat->dates->values[1] = 3; */
/*     dat->dna->list[0]->seq[0] = 'a'; */
/*     dat->dna->list[1]->seq[0] = 'a'; */
/*     dat->dna->list[2]->seq[0] = 't'; */
/*     printf("\n>>> Data <<<\n"); */
/*     print_data(dat); */


/*     /\* GET TIME SPAN *\/ */
/*     TIMESPAN = max_vec_int(dat->dates) - min_vec_int(dat->dates); */
/*     printf("\nTimespan is %d\n",TIMESPAN); */


/*     /\* CREATE AND INIT GENERATION TIME *\/ */
/*     gen = alloc_gentime(TIMESPAN, 5); */
/*     init_gentime(gen, 1, 1.0, 0.0, 0.0); */
/*     printf("\n>>> gentime info <<<\n"); */
/*     print_gentime(gen); */


/*      /\* CREATE AND INIT PARAMETERS *\/ */
/*     par = alloc_param(3); */
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


/*     /\* COMPUTE GENETIC DISTANCES *\/ */
/*     dnaInfo = compute_dna_distances(dat->dna); */
/*     printf("\n>>> DNA info <<<\n"); */
/*     print_dna_dist(dnaInfo); */


/*     /\* COMPUTE PRIORS *\/ */
/*     logPrior = logprior_all(par); */
/*     printf("\nPrior value (log): %.10f\n", logPrior); */

/*    /\* COMPUTE LIKELIHOOD *\/ */
/*     logLike = loglikelihood_all(dat, dnaInfo, gen, par); */
/*     printf("\nLog-likelihood value: %.10f\n", logLike); */

/*     /\* COMPUTE POSTERIOR *\/ */
/*     logPost = logposterior_all(dat, dnaInfo, gen, par); */
/*     printf("\nLog-posterior value: %.10f\n", logPost); */

/*     /\* /\\* RUNTIME TEST *\\/ *\/ */
/*     /\* int ITER=10e6, i; *\/ */
/*     /\* time_t t1, t2; *\/ */
/*     /\* time(&t1); *\/ */
/*     /\* printf("\nRuntime (%d computations of posterior): \n", ITER); *\/ */
/*     /\* for(i=0;i<ITER;i++){ *\/ */
/*     /\* 	logPost = logposterior_all(dat, dnaInfo, gen, par); *\/ */
/*     /\* } *\/ */
/*     /\* time(&t2); *\/ */
/*     /\* printf("\nellapsed time: %d seconds\n", (int) t2 - (int) t1); *\/ */


/*     /\* FREE / RETURN *\/ */
/*     free_data(dat); */
/*     free_gentime(gen); */
/*     free_dna_dist(dnaInfo); */
/*     free_param(par); */

/*     return 0; */
/* } */



/*
  gcc instructions

  gcc -o likelihood matvec.c genclasses.c structures.c init.c distances.c prior.c likelihood.c -lgsl -lgslcblas -Wall -g


  gcc -o likelihood matvec.c genclasses.c structures.c init.c distances.c prior.c likelihood.c -lgsl -lgslcblas -Wall -O3

 ./likelihood

  valgrind --leak-check=full -v likelihood

*/


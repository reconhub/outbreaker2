#include "common.h"
#include "structures.h"
#include "matvec.h"
#include "genclasses.h"
#include "distances.h"
#include "init.h"
#include "prior.h"
#include "likelihood.h"
#include "moves.h"




/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/
bool look_for_aba(param *par, data *dat){
  int B, A, ancesA;
  for(B=0;B<dat->n;B++){
    A = vec_int_i(par->alpha, B); /* A ancestor of B */
    if(A>-1){
      ancesA = vec_int_i(par->alpha, A); /*ancestor of A */
      if(ancesA==B){
	Rprintf("\nDetected pattern A->B->A: %d->%d->%d",A, B, A);
	Rprintf("\nCorrdesponding Tinf: A = %d,   B = %d", vec_int_i(par->Tinf,A), vec_int_i(par->Tinf,B));
	return TRUE;
      }
    }
  }
  return FALSE;
} /* end look_for_aba*/




/* SAMPLE KAPPA_I USING PROB */
/* sample a value of kappa_i according to respective proba */
/* of kappa_1, kappa_2, ..., kappa_maxK */
int choose_kappa_i(int T, gentime *gen, gsl_rng *rng){
    int i;
    double probVec[gen->maxK], sumKappa=0.0, cumSum=0.0, rnd;

    /* STANDARDIZE P(KAPPA_I) TO GET PROBA */
    for(i=1;i<gen->maxK;i++){
	probVec[i-1] = gentime_dens(gen, T, i);
	sumKappa += probVec[i-1];
    }

    for(i=0;i<gen->maxK;i++){
	probVec[i] = probVec[i]/sumKappa;
    }


    /* CHOOSE KAPPA_I WITH PROBA PROBVEC[I] */
    rnd = gsl_rng_uniform(rng);
    i=0;
    do{
	cumSum += probVec[i++];
    } while(rnd>=cumSum);

    if(i>gen->maxK){
      error("\n[in: moves.c->choose_kappa_i]\nInvalid value of kappa_i returned (%d, but maxK=%d).\n", i, gen->maxK);
      /* fprintf(stderr, "\n[in: moves.c->choose_kappa_i]\nInvalid value of kappa_i returned (%d, but maxK=%d).\n", i, gen->maxK); */
      /* exit(1); */
    }

    return i;
} /* end choose_kappa_i */





/*
   SAMPLE ALPHA_I USING PROB BASED ON MUTATIONS
   -> choose alpha from list of candidates
   -> 'alpha_i' is the most recent sampled ancestor of 'i'
   -> candidates = any earlier case
   -> returned value is the index of the proposed ancestor
   -> current version would allow for weights to be defined for candidates
     (but not used; note: removing this would make the procedure faster)
*/
int choose_alpha_i(int i, data *dat, param *par, mcmc_param *mcmcPar, gsl_rng *rng){
    int j, nCandidates, idOut, out;

    /* ESCAPE IF CURRENT ANCESTOR IS -1 */
    if(vec_int_i(par->alpha,i)<0) return -1;

    /* GET LIST OF CANDIDATES */
    nCandidates=0;

    /* candidates = infection time before T_i^inf*/
    for(j=0;j<dat->n;j++){
	if(vec_int_i(par->Tinf,j) < vec_int_i(par->Tinf,i)){
  	  /*   /\* debug *\/ */
	  /* if(j==i) Rprintf("\n### Messy stuff here: proposing %d->%d",j,i); */
	    /* store 'j' as a candidate */
	    mcmcPar->candid_ances->values[nCandidates] =  j;

	    /* store sampling weight for new candidate */
	    mcmcPar->candid_ances_proba->values[j] = 1.0;

	    /* increment number of candidates */
	    nCandidates++;
	} else {
	    mcmcPar->candid_ances_proba->values[j] = 0.0;
	}
    }


    /* RETURN PROPOSED ALPHA_I */
    /* no candidate = index case */
    if(nCandidates==0) return -1;

    /* one candidate */
    if(nCandidates==1) return vec_int_i(mcmcPar->candid_ances,0);

    /* >1 candidates */
    idOut = draw_multinom_censored(mcmcPar->candid_ances_proba, nCandidates, rng);
    out = vec_int_i(mcmcPar->candid_ances, idOut);

    return out;
} /* end choose_alpha_i */





/* find date of first imported case */
int find_date_first_import(data *dat, param *par){
  int i, out=0, counter=0;

  for(i=0;i<dat->n;i++){
    if(vec_int_i(par->alpha, i)<0){ /* imported case */
      if(counter==0){
	out = vec_int_i(par->Tinf, i);
      } else {
	if(out > vec_int_i(par->Tinf, i)) out = vec_int_i(par->Tinf, i);
      }
      counter++;
    }
  }
  return out;
}



/*
  =====
  MOVES
  =====
*/

/* MOVE VALUES OF MU1 */
void move_mu1(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, mcmc_param *mcmcPar, gsl_rng *rng){
    double logRatio=0.0, QCur, QTemp;

    /* GENERATE CANDIDATE VALUE FOR MU1 ACCORDING TO LOGNORMAL DISTRIBUTION */
    do{
	tempPar->mu1 = gsl_ran_lognormal(rng,log(currentPar->mu1),mcmcPar->sigma_mu1);
    } while(tempPar->mu1>1.0);
    /* which should be the same as:
    tempPar->mu1 = currentPar->mu1 * gsl_ran_lognormal(rng,0,mcmcPar->sigma_mu1); */

    /* other possibility for proposal */
    /* tempPar->mu1 += gsl_ran_gaussian(rng, mcmcPar->sigma_mu1); */
    /* if(tempPar->mu1 < 0.0) tempPar->mu1 = 0.0;*/

    /* ACCEPT / REJECT */
    /* only likelihood as priors are flat for mu1 */
    /* compute only genetic part as the epi part is unchanged */
    /* likelihoods */
    logRatio += loglikelihood_gen_all(dat, dnaInfo, tempPar, rng);
    logRatio -= loglikelihood_gen_all(dat, dnaInfo, currentPar, rng);
    /* priors */
    logRatio += logprior_mu1(tempPar);
    logRatio -= logprior_mu1(currentPar);

    /* add correction for MH truncated lognormal */
    QCur = gsl_cdf_gaussian_P(-log(currentPar->mu1),mcmcPar->sigma_mu1);
    QTemp = gsl_cdf_gaussian_P(-log(tempPar->mu1),mcmcPar->sigma_mu1);
    logRatio +=  log(tempPar->mu1) - log(currentPar->mu1); /* correction for lognormal */
    logRatio +=   log(QCur) - log(QTemp); /* correction for truncation (no values >1) */


    /* if p(new/old) > 1, accept new */
    if(logRatio>=0.0) {
	currentPar->mu1 = tempPar->mu1;
	mcmcPar->n_accept_mu1 += 1;
	/* printf("\nAccepting new value\n"); */
    } else { /* else accept new with proba (new/old) */
	if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
	    currentPar->mu1 = tempPar->mu1;
	    mcmcPar->n_accept_mu1 += 1;
	    /* printf("\nAccepting new value\n"); */
	} else { /* reject */
	    tempPar->mu1 = currentPar->mu1;
	    mcmcPar->n_reject_mu1 += 1;
	    /* printf("\nRejecting new value\n"); */
	}
    }

} /* end move_mu1 */





/* MOVE VALUES OF GAMMA */
void move_gamma(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, mcmc_param *mcmcPar, gsl_rng *rng){
    double logRatio=0.0;

    /* GENERATE CANDIDATE VALUE FOR GAMMA ACCORDING TO LOGNORMAL DISTRIBUTION */
    tempPar->gamma = gsl_ran_lognormal(rng,log(currentPar->gamma),mcmcPar->sigma_gamma);
    /* which should be the same as:
    tempPar->gamma = currentPar->gamma * gsl_ran_lognormal(rng,0,mcmcPar->sigma_gamma); */

      /* other possibility for proposal */
      /* tempPar->gamma += gsl_ran_gaussian(rng, mcmcPar->sigma_gamma);
    if(tempPar->gamma < 0.0) tempPar->gamma = 0.0;*/

    /* ACCEPT / REJECT */
    /* compute only genetic part as the epi part is unchanged */
    logRatio += loglikelihood_gen_all(dat, dnaInfo, tempPar, rng);
    logRatio -= loglikelihood_gen_all(dat, dnaInfo, currentPar, rng);

   /* add correction (MH) for lognormal proposal */
    logRatio += log(tempPar->gamma) - log(currentPar->gamma);

    /* compute the priors */
    logRatio += logprior_gamma(tempPar);
    logRatio -= logprior_gamma(currentPar);


    /* if p(new/old) > 1, accept new */
    if(logRatio>=0.0) {
	currentPar->gamma = tempPar->gamma;
	mcmcPar->n_accept_gamma += 1;
    } else { /* else accept new with proba (new/old) */
	if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
	    currentPar->gamma = tempPar->gamma;
	    mcmcPar->n_accept_gamma += 1;
	} else { /* reject */
	    tempPar->gamma = currentPar->gamma;
	    mcmcPar->n_reject_gamma += 1;
	}
    }

} /* end move_gamma*/





/* MOVE VALUES OF PI */
void move_pi(param *currentPar, param *tempPar, data *dat, mcmc_param *mcmcPar, gsl_rng *rng){
    double logRatio=0.0;
    double QCur, QTemp;

    /* GENERATE CANDIDATE VALUE FOR PI */
    /* HERE REPLACE WITH TRUNCATED LOGNORMAL (no values >1) )*/
    do
    {
    	tempPar->pi = gsl_ran_lognormal(rng,log(currentPar->pi),mcmcPar->sigma_pi);
    	/* which should be the same as: */
	/* tempPar->pi = currentPar->pi*gsl_ran_lognormal(rng,0,mcmcPar->sigma_pi); */
    } while(tempPar->pi>1.0);


    /* ACCEPT / REJECT */
    /* likelihood */
    logRatio += loglike_kappa_all(tempPar) - loglike_kappa_all(currentPar);

    /* prior */
    logRatio += logprior_pi(tempPar) - logprior_pi(currentPar);

    /* ADD CORRECTION FOR MH truncated lognormal */
    QCur = gsl_cdf_gaussian_P(-log(currentPar->pi),mcmcPar->sigma_pi);
    QTemp = gsl_cdf_gaussian_P(-log(tempPar->pi),mcmcPar->sigma_pi);
    logRatio +=  log(tempPar->pi) - log(currentPar->pi); /* correction for lognormal */
    logRatio +=   log(QCur) - log(QTemp); /* correction for truncation (no values >1) */

    /* if p(new/old) > 1, accept new */
    if(logRatio>=0.0) {
	currentPar->pi = tempPar->pi;
	mcmcPar->n_accept_pi += 1;
	/* printf("\nAccepting new value\n"); */
    } else { /* else accept new with proba (new/old) */
	if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
	    currentPar->pi = tempPar->pi;
	    mcmcPar->n_accept_pi += 1;
	    /* printf("\nAccepting new value\n"); */
	} else { /* reject */
	    tempPar->pi = currentPar->pi;
	    mcmcPar->n_reject_pi += 1;
	    /* printf("\nRejecting new value\n"); */
	}
    }

} /* end move_pi */






/* MOVE VALUES OF PHI */
void move_phi(param *currentPar, param *tempPar, data *dat, spatial_dist *spaInfo, mcmc_param *mcmcPar, gsl_rng *rng){
  double sumSameLoc=0.0, sumDiffLoc=0.0;
  int i, ances;

  /* USE GIBBS SAMPLER */
  for(i=0;i < dat->n;i++){
    ances = vec_int_i(currentPar->alpha,i);

    if(ances>=0){
      if(vec_int_i(dat->locations,i)==vec_int_i(dat->locations,ances)){
	sumSameLoc += 1.0;
      } else {
	sumDiffLoc += 1.0;
      }
    }
  }

  tempPar->phi = gsl_ran_beta(rng, currentPar->phi_param1 + sumSameLoc, currentPar->phi_param2 + sumDiffLoc);
  currentPar->phi = tempPar->phi;

    /* /\* GENERATE CANDIDATE VALUE FOR PHI *\/ */
    /* /\* HERE REPLACE WITH TRUNCATED LOGNORMAL (no values >1) )*\/ */
    /* do */
    /* { */
    /* 	tempPar->phi = gsl_ran_lognormal(rng,log(currentPar->phi),mcmcPar->sigma_phi); */
    /* 	/\* which should be the same as: *\/ */
    /* 	/\* tempPar->phi = currentPar->phi*gsl_ran_lognormal(rng,0,mcmcPar->sigma_phi); *\/ */
    /* } while(tempPar->phi>1.0); */


    /* /\* ACCEPT / REJECT *\/ */
    /* /\* likelihood *\/ */
    /* 	logRatio += loglikelihood_spa_all(dat, spaInfo, tempPar, rng) - loglikelihood_spa_all(dat, spaInfo, currentPar, rng); */

    /* /\* prior *\/ */
    /* logRatio += logprior_phi(tempPar) - logprior_phi(currentPar); */

    /* /\* ADD CORRECTION FOR MH truncated lognormal *\/ */
    /* QCur = gsl_cdf_gaussian_P(-log(currentPar->phi),mcmcPar->sigma_phi); */
    /* QTemp = gsl_cdf_gaussian_P(-log(tempPar->phi),mcmcPar->sigma_phi); */
    /* logRatio +=  log(tempPar->phi) - log(currentPar->phi); /\* correction for lognormal *\/ */
    /* logRatio +=   log(QCur) - log(QTemp); /\* correction for truncation (no values >1) *\/ */

    /* /\* if p(new/old) > 1, accept new *\/ */
    /* if(logRatio>=0.0) { */
    /* 	currentPar->phi = tempPar->phi; */
    /* 	mcmcPar->n_accept_phi += 1; */
    /* 	/\* printf("\nAccepting new value\n"); *\/ */
    /* } else { /\* else accept new with proba (new/old) *\/ */
    /* 	if(log(gsl_rng_uniform(rng)) <= logRatio){ /\* accept *\/ */
    /* 	    currentPar->phi = tempPar->phi; */
    /* 	    mcmcPar->n_accept_phi += 1; */
    /* 	    /\* printf("\nAccepting new value\n"); *\/ */
    /* 	} else { /\* reject *\/ */
    /* 	    tempPar->phi = currentPar->phi; */
    /* 	    mcmcPar->n_reject_phi += 1; */
    /* 	    /\* printf("\nRejecting new value\n"); *\/ */
    /* 	} */
    /* } */

} /* end move_phi */






/* MOVE VALUES OF SPA_PARAM1 */
void move_spa1(param *currentPar, param *tempPar, data *dat, spatial_dist *spaInfo, mcmc_param *mcmcPar, gsl_rng *rng){
  double logRatio=0.0;

  /* MOVEMENTS AND CORRECTIONS ARE MODEL-DEPENDENT*/
  /* spatial model 0: no spatial info */
  if(currentPar->spa_model == 0) return;

  /* spatial model 1: exponential diffusion */
  if(currentPar->spa_model == 1){
    /* propose new value */

    tempPar->spa_param1 = gsl_ran_lognormal(rng, log(currentPar->spa_param1), mcmcPar->sigma_spa1);
    /* ACCEPT / REJECT */
    /* compute only spatial part of likelihood as the rest is unchanged */
    logRatio += loglikelihood_spa_all(dat, spaInfo, tempPar, rng);
    logRatio -= loglikelihood_spa_all(dat, spaInfo, currentPar, rng);

    /* add correction (MH) for lognormal proposal if needed */
    if(currentPar->spa_model>0){
      logRatio += log(tempPar->spa_param1) - log(currentPar->spa_param1);
    }

    /* compute the priors */
    logRatio += logprior_spa1(tempPar);
    logRatio -= logprior_spa1(currentPar);

    /* if p(new/old) > 1, accept new */
    if(logRatio>=0.0) {
      currentPar->spa_param1 = tempPar->spa_param1;
      mcmcPar->n_accept_spa1 += 1;
    } else { /* else accept new with proba (new/old) */
      if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
	currentPar->spa_param1 = tempPar->spa_param1;
	mcmcPar->n_accept_spa1 += 1;
      } else { /* reject */
	tempPar->spa_param1 = currentPar->spa_param1;
	mcmcPar->n_reject_spa1 += 1;
      }
    }
  }
} /* end move_spa1 */




/* MOVE INFECTION DATES (T_inf) */
void move_Tinf(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng){
    double logRatio=0.0;
    int i, toMove = 0;

    /* DETERMINE WHICH Tinf_i TO MOVE */
    sample_vec_int(mcmcPar->all_idx, mcmcPar->idx_move_Tinf, FALSE, rng);

    /* MOVE EACH Tinf_i IN TURN */
    for(i=0;i<mcmcPar->idx_move_Tinf->length;i++){
	toMove = vec_int_i(mcmcPar->idx_move_Tinf,i);

	/* move i-th Tinf */
	tempPar->Tinf->values[toMove] += (gsl_rng_uniform(rng) >= 0.5 ? 1 : -1) * gsl_ran_poisson(rng, 1);

	/* MAY NEED TO CHANGE THIS AND ADD CORRECTION */
	/* constraint: Tinf_i < t_i */
	if(!vec_int_i(tempPar->Tinf,toMove) < vec_int_i(dat->dates,toMove)) tempPar->Tinf->values[toMove] = vec_int_i(dat->dates,toMove)-1;
	/* constraint: Tinf_i >= -truncW */
	if(vec_int_i(tempPar->Tinf,toMove) < -gen->truncW) tempPar->Tinf->values[toMove] = -gen->truncW;

	/* PROCEED TO ACCEPT/REJECT ONLY IF TINF HAS CHANGED */
	if(vec_int_i(tempPar->Tinf,toMove) != vec_int_i(currentPar->Tinf,toMove)){
	    /* ACCEPT/REJECT STEP */
	    /* compute the likelihood (no priors for Tinf) */
	    /* compute only local changes in the likelihood */
	    logRatio = loglikelihood_local_i(toMove, dat, dnaInfo, spaInfo, gen, tempPar, rng) - loglikelihood_local_i(toMove, dat, dnaInfo, spaInfo, gen, currentPar, rng);
	    /* logRatio = loglikelihood_all(dat, dnaInfo, spaInfo, gen, tempPar, rng) - loglikelihood_all(dat, dnaInfo, spaInfo, gen, currentPar, rng); */


	    /* if p(new/old) > 1, accept new */
	    if(logRatio>=0.0) {
		/* printf("\nTinf_%d: accepting automatically move from %d to %d (respective loglike:%f and %f)\n",toMove+1, vec_int_i(currentPar->Tinf,toMove), vec_int_i(tempPar->Tinf,toMove), ll1, ll2); */
		/* fflush(stdout); */

		currentPar->Tinf->values[toMove] = vec_int_i(tempPar->Tinf,toMove);
		mcmcPar->n_accept_Tinf += 1;
	    } else { /* else accept new with proba (new/old) */
		if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
		/*     printf("\nTinf_%d: accepting move from %d to %d (respective loglike:%f and %f)\n",toMove+1, vec_int_i(currentPar->Tinf,toMove), vec_int_i(tempPar->Tinf,toMove), ll1, ll2); */
		/* fflush(stdout); */

		    currentPar->Tinf->values[toMove] = vec_int_i(tempPar->Tinf,toMove);
		    mcmcPar->n_accept_Tinf += 1;
		} else { /* reject */
		    tempPar->Tinf->values[toMove] = vec_int_i(currentPar->Tinf,toMove);
		    mcmcPar->n_reject_Tinf += 1;
		}
	    }
	} /* end if Tinf has changed */
    } /* end for each indiv to move */
} /* end move_Tinf*/







/* MOVE INFECTION DATES, NB OF GENERATIONS, AND ANCESTRIES */
void move_Tinf_alpha_kappa(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng){
  int i, j, toMove=0, nbCandidCurrent=0, nbCandidTemp=0, nbDaysCurrent=0, nbDaysTemp=0, firstImported=0, ances;
  double logRatio = 0.0, correcRatio = 0.0, ll1 = 0.0, ll2 = 0.0;

  /* DEBUGGING */
  /* bool checkLL = check_loglikelihood_all(dat, dnaInfo, spaInfo, gen, currentPar, rng); */

  /* if(!checkLL) Rprintf("\n-Inf LL detected when entering move_Tinf_alpha_kappa"); */

  /* bool checkABA = look_for_aba(currentPar, dat); */
  /* bool oldCheckABA = checkABA; */
  /* if(checkABA) Rprintf("\nABA detected when entering move_Tinf_alpha_kappa"); */

  /* DETERMINE WHICH INDIVIDUAL TO MOVE */
  draw_vec_int_multinom(mcmcPar->all_idx, mcmcPar->idx_move_alpha, mcmcPar->move_alpha, rng);

  /* MOVE Tinf, kappa, alpha FOR EACH CHOSEN INDIVIDUAL */
  for(i=0;i<mcmcPar->idx_move_alpha->length;i++){
    /* initialize correction for asymetric proposal */
    correcRatio=0.0;

    /* CURRENT INDIVIDUAL TO MOVE */
    toMove = vec_int_i(mcmcPar->idx_move_alpha,i);

    /* make sure that we don't move imported cases */
    ances = vec_int_i(currentPar->alpha, toMove);
    if(ances>-1){

      /* MOVE Tinf UNLESS USER DISABLED THIS MOVE */
      if(mcmcPar->move_Tinf){
	/* find first imported case */
	firstImported = find_date_first_import(dat, currentPar);

	/* move i-th Tinf */
	tempPar->Tinf->values[toMove] += (gsl_rng_uniform(rng) >= 0.5 ? 1 : -1) * gsl_ran_poisson(rng, 1);

	/* constraint: Tinf_i < t_i */
	if(vec_int_i(tempPar->Tinf,toMove) >= vec_int_i(dat->dates,toMove)) tempPar->Tinf->values[toMove] = vec_int_i(dat->dates,toMove)-1;

	/* constraint: Tinf_i > first imported */
	if(vec_int_i(tempPar->Tinf,toMove) <= firstImported) tempPar->Tinf->values[toMove] = firstImported+1;
      }


      /* MOVE ALPHA IF MOVEABLE  */
      if(vec_double_i(mcmcPar->move_alpha,toMove)>0.0){
	/* check number of possible ancestors before/after move */
	nbCandidCurrent=0;
	nbCandidTemp=0;
	for(j=0;j<dat->n;j++){
	  if(vec_int_i(currentPar->Tinf,j) < vec_int_i(currentPar->Tinf,toMove))
	    nbCandidCurrent++;
	  if(vec_int_i(tempPar->Tinf,j) < vec_int_i(tempPar->Tinf,toMove))
	    nbCandidTemp++;
	}

	/* propose new alpha */
	tempPar->alpha->values[toMove] = choose_alpha_i(toMove, dat, tempPar, mcmcPar, rng);

	/* compute correction factor */
	/* log[(1/nbCandid) / (1/nbCandid*)] */
	/* = log(1) - log(nbCandid) - log(1) + log(= nbCandid*) */
	/* = log(nbCandid*) - log(nbCandid) */
	/* if(nbCandidCurrent<1) { */
	/*   warning("\n[in move_Tinf_alpha_kappa]: nbCandidCurrent is zero and will result in a dodgy move; are you trying to move an imported case?"); */
	/*   Rprintf("\nnbCandidCurrent is zero"); */
	/*   Rprintf("\ndate of first imported case is %d",firstImported); */
	/*   Rprintf("\nproposed Tinf for %d is %d", toMove, tempPar->Tinf->values[toMove]); */
	/*   Rprintf("\nThis happened while moving %d:",toMove); */

	/*   Rprintf("\ntrying to move from %d->%d to %d->%d (respective loglike: old:%f  new%f   correc:%f)\n",vec_int_i(currentPar->alpha,toMove), toMove, vec_int_i(tempPar->alpha,toMove), toMove, ll1, ll2, correcRatio); */

	/*   Rprintf("\n====== Current param where: ======="); */
	/*   print_param(currentPar); */

	/*   Rprintf("\n\n====== New param are: ======"); */
	/*   print_param(tempPar); */
	/*   getchar(); */
	/* } */
	correcRatio += log(nbCandidTemp) - log(nbCandidCurrent);
      }


      /* MOVE KAPPA */
      /* if not imported and moveable */
      if(vec_int_i(tempPar->alpha,toMove)>=0 && vec_double_i(mcmcPar->move_kappa,toMove)>0.0){
	/* propose new, 'intelligent' kappa */
	nbDaysTemp = vec_int_i(tempPar->Tinf,toMove) - vec_int_i(tempPar->Tinf, vec_int_i(tempPar->alpha,toMove));
	tempPar->kappa->values[toMove] = choose_kappa_i(nbDaysTemp, gen, rng);

	/* compute correction factor */
	/* formula: */
	/* log[w^(kappa_i)(T_i^inf - T_{alpha_i}^inf) / w^(kappa_i*)(T_i^inf* - T_{alpha_i*}^inf)] */
	/* = log(w^(kappa_i)(T_i^inf - T_{alpha_i}^inf) - log(w^(kappa_i*)(T_i^inf* - T_{alpha_i*}^inf)) */
	nbDaysCurrent = vec_int_i(currentPar->Tinf,toMove) - vec_int_i(currentPar->Tinf, vec_int_i(currentPar->alpha,toMove));
	correcRatio += log(gentime_dens(gen, nbDaysCurrent, vec_int_i(currentPar->kappa,toMove)));
	correcRatio -= log(gentime_dens(gen, nbDaysTemp, vec_int_i(tempPar->kappa,toMove)));
      }


      /* ACCEPT/REJECT STEP */
      /* compute the likelihood ratio */
      /* ll1 = loglikelihood_all(dat, dnaInfo, spaInfo, gen, currentPar, rng); */
      /* ll2 = loglikelihood_all(dat, dnaInfo, spaInfo, gen, tempPar, rng); */
      /* compute only local changes in the likelihood */
      ll1 = loglikelihood_local_i(toMove, dat, dnaInfo, spaInfo, gen, currentPar, rng);
      ll2 = loglikelihood_local_i(toMove, dat, dnaInfo, spaInfo, gen, tempPar, rng);
      logRatio = ll2 - ll1;
      filter_logprob(&correcRatio);
      logRatio += correcRatio;

      /* /\* MH correction *\/ */
      /* /\* like ratio x ( Pmove(current)/Pmove(temp) ) *\/ */


      /* if p(new/old) > 1, accept new */
      if(logRatio>=0.0) {
	/* /\* DEBUGGING *\/ */
	/* if(check_loglikelihood_all(dat, dnaInfo, spaInfo, gen, currentPar, rng) && !check_loglikelihood_all(dat, dnaInfo, spaInfo, gen, tempPar, rng)) { */
	/*   Rprintf("\nMovement to -Inf LL during move_Tinf_alpha_kappa"); */
	/*   Rprintf("\nThis happened while moving %d",toMove); */

	/*   Rprintf("\naccepting automatically move from %d->%d to %d->%d (respective loglike: old:%f  new%f   correc:%f)\n",vec_int_i(currentPar->alpha,toMove), toMove, vec_int_i(tempPar->alpha,toMove), toMove, ll1, ll2, correcRatio); */

	/*   Rprintf("\n====== Current param where: ======="); */
	/*   print_param(currentPar); */

	/*   Rprintf("\n\n====== New param are: ======"); */
	/*   print_param(tempPar); */
	/*   getchar(); */
	/* } */

	/* /\* DEBUGGING *\/ */
	/* if(!look_for_aba(currentPar, dat) && look_for_aba(tempPar, dat)) { */
	/*   Rprintf("\nABA created during move_Tinf_alpha_kappa"); */
	/*   Rprintf("\nThis happened while moving %d",toMove); */

	/*   Rprintf("\naccepting automatically move from %d->%d to %d->%d (respective loglike: old:%f  new%f   correc:%f)\n",vec_int_i(currentPar->alpha,toMove), toMove, vec_int_i(tempPar->alpha,toMove), toMove, ll1, ll2, correcRatio); */

	/*   Rprintf("\n====== Current param where: ======="); */
	/*   print_param(currentPar); */

	/*   Rprintf("\n\n====== New param are: ======"); */
	/*   print_param(tempPar); */
	/*   getchar(); */

	/* } */

	currentPar->Tinf->values[toMove] = vec_int_i(tempPar->Tinf,toMove);
	currentPar->alpha->values[toMove] = vec_int_i(tempPar->alpha,toMove);
	currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove);
	mcmcPar->n_accept_Tinf += 1;
	mcmcPar->n_accept_alpha += 1;
	mcmcPar->n_accept_kappa += 1;
      } else { /* else accept new with proba (new/old) */
	if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
	  /* /\* DEBUGGING *\/ */
	  /* if(check_loglikelihood_all(dat, dnaInfo, spaInfo, gen, currentPar, rng) && !check_loglikelihood_all(dat, dnaInfo, spaInfo, gen, tempPar, rng)) { */
	  /*   Rprintf("\nMovement to -Inf LL during move_Tinf_alpha_kappa"); */
	  /*   Rprintf("\nThis happened while moving %d",toMove); */

	  /*   Rprintf("\naccepting automatically move from %d->%d to %d->%d (respective loglike: old:%f  new%f   correc:%f)\n",vec_int_i(currentPar->alpha,toMove), toMove, vec_int_i(tempPar->alpha,toMove), toMove, ll1, ll2, correcRatio); */

	  /*   Rprintf("\n====== Current param where: ======="); */
	  /*   print_param(currentPar); */

	  /*   Rprintf("\n\n====== New param are: ======"); */
	  /*   print_param(tempPar); */
	  /*   getchar(); */
	  /* } */

	  /* if(!look_for_aba(currentPar, dat) && look_for_aba(tempPar, dat)) { */
	  /*   Rprintf("\nABA created during move_Tinf_alpha_kappa"); */
	  /*   Rprintf("\nThis happened while moving %d",toMove); */

	  /*   Rprintf("\naccepting automatically move from %d->%d to %d->%d (respective loglike:%f and %f)\n",vec_int_i(currentPar->alpha,toMove), toMove, vec_int_i(tempPar->alpha,toMove), toMove, ll1, ll2); */

	  /*   Rprintf("\nCurrent param where:"); */
	  /*   print_param(currentPar); */

	  /*   Rprintf("\nNew param are:"); */
	  /*   print_param(tempPar); */
	  /*   getchar(); */
	  /* } */

	  /* /\* debugging *\/ */
	  /* printf("\naccepting move from %d->%d to %d->%d (respective loglike:%f and %f)\n",vec_int_i(currentPar->alpha,toMove), toMove+1, vec_int_i(tempPar->alpha,toMove), toMove+1, ll1, ll2); */
	  /* fflush(stdout); */
	  currentPar->Tinf->values[toMove] = vec_int_i(tempPar->Tinf,toMove);
	  currentPar->alpha->values[toMove] = vec_int_i(tempPar->alpha,toMove);
	  currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove);
	  mcmcPar->n_accept_Tinf += 1;
	  mcmcPar->n_accept_alpha += 1;
	  mcmcPar->n_accept_kappa += 1;
	} else { /* reject */
	  tempPar->Tinf->values[toMove] = vec_int_i(currentPar->Tinf,toMove);
	  tempPar->alpha->values[toMove] = vec_int_i(currentPar->alpha,toMove);
	  tempPar->kappa->values[toMove] = vec_int_i(currentPar->kappa,toMove);
	  mcmcPar->n_reject_Tinf += 1;
	  mcmcPar->n_reject_alpha += 1;
	}
      } /* end  ACCEPT/REJECT STEP */
    } /* end if ances > -1 */
  } /* end for loop (for all 'i' to move) */

    /* /\* DEBUGGING *\/ */
    /* checkABA = look_for_aba(currentPar, dat); */
    /* if(!oldCheckABA && checkABA) { */
    /*   Rprintf("\nABA created during move_Tinf_alpha_kappa"); */
    /*   getchar(); */
    /* } else if(checkABA) Rprintf("\nABA detected when leaving move_Tinf_alpha_kappa"); */


} /* end move_Tinf_alpha_kappa */








/* swap ancestries in the chain x->A->B */
/* B is picked randomly, then A and B are swapped so that: */
/* x->A->B becomes x->B->A */
/* all descendents of B become descendents of A */
/* all descendents of A become descendents of B */
void swap_ancestries(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng){
  int i, j, x, A, B;
  double logRatio = 0.0;


  /* /\* DEBUGGING *\/ */
  /*   bool checkABA = look_for_aba(currentPar, dat); */
  /*   bool oldCheckABA = checkABA; */
  /*   if(checkABA) Rprintf("\nABA detected when entering swap_ancestries"); */

  /* DETERMINE WHICH INDIVIDUAL TO MOVE */
  draw_vec_int_multinom(mcmcPar->all_idx, mcmcPar->idx_move_alpha, mcmcPar->move_alpha, rng);

  /* MOVE Tinf, kappa, alpha FOR EACH CHOSEN INDIVIDUAL */

  /* Rprintf("\n Forcing Tinf to [-2;0;1]"); */
  /* currentPar->Tinf->values[0]=-2; */
  /* currentPar->Tinf->values[1]=0; */
  /* currentPar->Tinf->values[2]=1; */

  for(i=0;i<mcmcPar->idx_move_alpha->length;i++){

    /* CURRENT INDIVIDUAL TO MOVE */
    /* reminder: x->A->B; x, A, B are indices of cases */
    /* 'B': case whose ancestry is inverted */
    /* 'A': vec_int_i(tempPar->alpha,B) */
    /* 'x': vec_int_i(tempPar->alpha, A) */
    B = vec_int_i(mcmcPar->idx_move_alpha, i);
    A = vec_int_i(currentPar->alpha, B);
 
    /* proceed only if B isn't imported */
    if(A>-1){
      x = vec_int_i(currentPar->alpha, A);

      /* SWAP ONLY IF: A AND B CAN MOVE, AND A ISN'T IMPORTED */
      /* i.e., A>-1 && x>-1 */
      if(vec_double_i(mcmcPar->move_alpha, A)>0.0 && vec_double_i(mcmcPar->move_alpha, B)>0.0 && x>-1){

	/* SWAP ANCESTRIES */
	tempPar->alpha->values[A] = B; /* (x->A) changes to (B->A) */
	tempPar->alpha->values[B] = x; /* (A->B) changes to (x->B) */

	/* SWAP TINF */
	tempPar->Tinf->values[A] = vec_int_i(currentPar->Tinf,B);
	tempPar->Tinf->values[B] = vec_int_i(currentPar->Tinf,A);

	/* ALL DESCENDENTS OF B BECOME DESCENDENTS OF A */
	/* ALL DESCENDENTS OF A BECOME DESCENDENTS OF B */
	for(j=0;j<dat->n;j++){
	  /* if case 'j' can move... */
	  if(vec_double_i(mcmcPar->move_alpha, j)>0.0){
	    /* if this case was descendent of B, it becomes descendent of A */
	    if(j!=A && vec_int_i(currentPar->alpha, j)==B){
	      tempPar->alpha->values[j] = A;
	      /* /\* DEBUGGING *\/ */
	      /* if(j==A) { */
	      /*   Rprintf("\nHere's a mess: j==A"); */
	      /*   Rprintf("\nThis happened while... "); */
	      /*   Rprintf("\nswapping %d->%d->%d to %d->%d->%d \n", x,A,B,x,B,A); */
	      /*   Rprintf("\nTinf vector before:  "); */
	      /*   print_vec_int(currentPar->Tinf); */
	      /*   Rprintf("\nTinf vector proposed:  "); */
	      /*   print_vec_int(tempPar->Tinf); */

	      /*   Rprintf("\nAlpha vector before:  "); */
	      /*   print_vec_int(currentPar->alpha); */
	      /*   Rprintf("\nAlpha vector proposed:  "); */
	      /*   print_vec_int(tempPar->alpha); */

	      /*   getchar(); */
	      /* } */
	      /* ...and descendents of A becomes descendent of B (except for B!!) */
	    } else if(j!=B && vec_int_i(currentPar->alpha, j)==A){
	      tempPar->alpha->values[j] = B;
	      /* /\* DEBUGGING *\/ */
	      /* if(j==B) Rprintf("Here's the mess: j==B"); */
	    }
	  }
	}


	/* ACCEPT/REJECT STEP */
	/* compute the likelihood ratio */
	/* logRatio = loglikelihood_all(dat, dnaInfo, spaInfo, gen, tempPar, rng) - loglikelihood_all(dat, dnaInfo, spaInfo, gen, currentPar, rng); */
	/* need to be careful with local likelihood: */
	/* in 'current' config, LL(B) is counted twice as because of A->B */
	/* in 'temp' config, LL(A) is counted twice as because of B->A */
	logRatio = loglikelihood_local_i(A, dat, dnaInfo, spaInfo, gen, tempPar, rng) +
	  loglikelihood_local_i(B, dat, dnaInfo, spaInfo, gen, tempPar, rng) -
	  loglikelihood_local_i(A, dat, dnaInfo, spaInfo, gen, currentPar, rng) -
	  loglikelihood_local_i(B, dat, dnaInfo, spaInfo, gen, currentPar, rng) -
	  loglikelihood_i(A, dat, dnaInfo, spaInfo, gen, tempPar, rng) +
	  loglikelihood_i(B, dat, dnaInfo, spaInfo, gen, currentPar, rng);
	/* ll2 = loglikelihood_all(dat, dnaInfo, spaInfo, gen, tempPar, rng); */
	/* ll1 = loglikelihood_all(dat, dnaInfo, spaInfo, gen, currentPar, rng); */
	/* logRatio = ll2 - ll1; */


	/* if(!check_loglikelihood_all(dat, dnaInfo, spaInfo, gen, tempPar, rng)){ */
	/* /\* Rprintf("\nThis happened while... "); *\/ */
	/* /\* Rprintf("\nswapping %d->%d->%d to %d->%d->%d \n", x,A,B,x,B,A); *\/ */
	/* /\* Rprintf("\nTinf vector before:  "); *\/ */
	/* /\* print_vec_int(currentPar->Tinf); *\/ */
	/* /\* Rprintf("\nTinf vector proposed:  "); *\/ */
	/* /\* print_vec_int(tempPar->Tinf); *\/ */

	/* /\* Rprintf("\nAlpha vector before:  "); *\/ */
	/* /\* print_vec_int(currentPar->alpha); *\/ */
	/* /\* Rprintf("\nAlpha vector proposed:  "); *\/ */
	/* /\* print_vec_int(tempPar->alpha); *\/ */

	/* } */
	/* Rprintf("\n"); */
	/* Rprintf("(new LL: %f) (old LL %f)", ll2, ll1); */

	/* if p(new/old) > 1, accept new */
	if(logRatio>=0.0) {
	  /* Rprintf("x"); */
	  copy_param(tempPar, currentPar);
	  mcmcPar->n_accept_Tinf += 1;
	  mcmcPar->n_accept_alpha += 1;
	} else { /* else accept new with proba (new/old) */
	  if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
	    /* Rprintf("x"); */
	    /* Rprintf("(new LL: %f) (old LL %f)", ll2, ll1); */
	    /* fflush(stdout); */
	    /* Rprintf("...accepted, by chance");fflush(stdout); */

	    /* /\* debugging *\/ */
	    /* printf("\naccepting move from %d->%d to %d->%d (respective loglike:%f and %f)\n",vec_int_i(currentPar->alpha,toMove), toMove+1, vec_int_i(tempPar->alpha,toMove), toMove+1, ll1, ll2); */
	    /* fflush(stdout); */

	    copy_param(tempPar, currentPar);
	    mcmcPar->n_accept_Tinf += 1;
	    mcmcPar->n_accept_alpha += 1;
	  } else { /* reject */
	    /* Rprintf("."); */

	    copy_param(currentPar, tempPar);
	    mcmcPar->n_reject_Tinf += 1;
	    mcmcPar->n_reject_alpha += 1;
	  }
	} /* end  ACCEPT/REJECT STEP */
      }
    }
    /* Rprintf("\n"); */
  } /* end for loop for all moved individuals */

  /* /\* DEBUGGING *\/ */
  /*   checkABA = look_for_aba(currentPar, dat); */
  /*   if(!oldCheckABA && checkABA) { */
  /*     Rprintf("\nABA created during swap_ancestries"); */
  /*     getchar(); */
  /*   } else if(checkABA) Rprintf("\nABA detected when leaving swap_ancestries"); */

} /* end swap_ancestries */

















/* NO LONGER USED */


/* /\* MOVE VALUES OF ALPHA AND KAPPA *\/ */
/* void move_alpha_kappa(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng){ */
/*     int i, toMove=0, temp; */
/*     double logRatio = 0.0; */


/*     /\* DETERMINE WHICH ALPHA_I/KAPPA_I TO MOVE *\/ */
/*     /\* sample_vec_int(mcmcPar->all_idx, mcmcPar->idx_move_alpha, FALSE, rng); *\/ */
/*     draw_vec_int_multinom(mcmcPar->all_idx, mcmcPar->idx_move_alpha, mcmcPar->move_alpha, rng); */

/*     /\* MOVE EACH ALPHA_I/KAPPA_I TO MOVE *\/ */
/*     /\* need to propose a new kappa at the same time *\/ */
/*     for(i=0;i<mcmcPar->idx_move_alpha->length;i++){ */
/* 	toMove = vec_int_i(mcmcPar->idx_move_alpha,i); */

/* 	/\* MOVE KAPPA *\/ */
/* 	/\* if not imported and moveable *\/ */
/* 	if(vec_int_i(tempPar->alpha,toMove)>=0 && vec_double_i(mcmcPar->move_kappa,toMove)>0.0){ */
/* 	    /\* movement *\/ */
/* 	    temp = tempPar->kappa->values[toMove] + (gsl_rng_uniform(rng) >= 0.5 ? 1 : -1); */

/* 	    /\* needs to be on [1;maxK]*\/ */
/* 	    if(temp < 1) { */
/* 		temp = 1; */
/* 	    } else if(temp>gen->maxK){ */
/* 		temp = gen->maxK; */
/* 	    } */
/* 	    /\* store new value of kappa_i *\/ */
/* 	    tempPar->kappa->values[toMove] = temp; */
/* 	} */

/* 	/\* MOVE ALPHA *\/ */
/* 	/\* if moveable *\/ */
/* 	if(vec_double_i(mcmcPar->move_alpha,toMove)>0.0){ */
/* 	    tempPar->alpha->values[toMove] = choose_alpha_i(toMove, dat, currentPar, mcmcPar, rng); */
/* 	} */

/* 	/\* ACCEPT/REJECT STEP *\/ */
/* 	/\* compute the likelihood *\/ */
/* 	/\* logRatio = loglikelihood_all(dat, dnaInfo, gen, tempPar, rng) - loglikelihood_all(dat, dnaInfo, gen, currentPar, rng); *\/ */
/* 	logRatio = loglikelihood_i(toMove, dat, dnaInfo, spaInfo, gen, tempPar, rng) - loglikelihood_i(toMove, dat, dnaInfo, spaInfo, gen, currentPar, rng); */

/* 	/\* /\\* MH correction *\\/ *\/ */
/* 	/\* /\\* like ratio x ( Pmove(current)/Pmove(temp) ) *\\/ *\/ */

/* 	/\* /\\* debugging *\\/ *\/ */
/* 	/\* ll1=loglikelihood_all(dat, dnaInfo, gen, currentPar, rng); *\/ */
/* 	/\* ll2=loglikelihood_all(dat, dnaInfo, gen, tempPar, rng); *\/ */

/* 	/\* if p(new/old) > 1, accept new *\/ */
/* 	if(logRatio>=0.0) { */
/* 	    /\* /\\* debugging *\\/ *\/ */
/* 	    /\* printf("\naccepting automatically move from %d->%d to %d->%d (respective loglike:%f and %f)\n",vec_int_i(currentPar->alpha,toMove), toMove+1, vec_int_i(tempPar->alpha,toMove), toMove+1, ll1, ll2); *\/ */
/* 	    /\* fflush(stdout); *\/ */

/* 	    currentPar->alpha->values[toMove] = vec_int_i(tempPar->alpha,toMove); */
/* 	    currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove); */
/* 	    /\* mcmcPar->Pmove_alpha_old = mcmcPar->Pmove_alpha_new; *\/ */
/* 	    mcmcPar->n_accept_alpha += 1; */
/* 	} else { /\* else accept new with proba (new/old) *\/ */
/* 	    if(log(gsl_rng_uniform(rng)) <= logRatio){ /\* accept *\/ */
/* 		/\* /\\* debugging *\\/ *\/ */
/* 		/\* printf("\naccepting move from %d->%d to %d->%d (respective loglike:%f and %f)\n",vec_int_i(currentPar->alpha,toMove), toMove+1, vec_int_i(tempPar->alpha,toMove), toMove+1, ll1, ll2); *\/ */
/* 		/\* fflush(stdout); *\/ */

/* 		currentPar->alpha->values[toMove] = vec_int_i(tempPar->alpha,toMove); */
/* 		currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove); */
/* 		/\* mcmcPar->Pmove_alpha_old = mcmcPar->Pmove_alpha_new; *\/ */
/* 		mcmcPar->n_accept_alpha += 1; */
/* 	    } else { /\* reject *\/ */
/* 		tempPar->alpha->values[toMove] = vec_int_i(currentPar->alpha,toMove); */
/* 		tempPar->kappa->values[toMove] = vec_int_i(currentPar->kappa,toMove); */
/* 		mcmcPar->n_reject_alpha += 1; */
/* 	    } */
/* 	} /\* end  ACCEPT/REJECT STEP *\/ */
/*     } /\* end for loop (for all 'i' to move) *\/ */
/* } /\* end move_alpha_kappa *\/ */





/* /\* MOVE VALUES OF ALPHA *\/ */
/* void move_alpha(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng){ */
/*     int i, toMove=0, T; */
/*     double logRatio = 0.0; */


/*     /\* DETERMINE WHICH ALPHA_I TO MOVE *\/ */
/*     /\* sample_vec_int(mcmcPar->all_idx, mcmcPar->idx_move_alpha, FALSE, rng); *\/ */
/*     draw_vec_int_multinom(mcmcPar->all_idx, mcmcPar->idx_move_alpha, mcmcPar->move_alpha, rng); */

/*     /\* MOVE EACH ALPHA_I TO MOVE *\/ */
/*     /\* need to propose a new kappa at the same time *\/ */
/*     for(i=0;i<mcmcPar->idx_move_alpha->length;i++){ */
/* 	toMove = vec_int_i(mcmcPar->idx_move_alpha,i); */

/* 	/\* MOVE ALPHA *\/ */
/* 	tempPar->alpha->values[toMove] = choose_alpha_i(toMove, dat, currentPar, mcmcPar, rng); */

/* 	/\* MOVE KAPPA *\/ */
/* 	/\* proceed only if alpha has changed *\/ */
/* 	if(vec_int_i(tempPar->alpha,toMove) != vec_int_i(currentPar->alpha,toMove)){ */

/* 	    if(vec_int_i(tempPar->alpha,toMove)>=0){ */
/* 		/\* move kappa_i - multinomial *\/ */
/* 		/\* T: Tinf_i - Tinf_ances *\/ */
/* 		T = vec_int_i(tempPar->Tinf,toMove) - vec_int_i(tempPar->Tinf, tempPar->alpha->values[toMove]); */
/* 		tempPar->kappa->values[toMove] = choose_kappa_i(T, gen, rng); */
/* 	    } else { /\* default value for imported cases *\/ */
/* 		tempPar->kappa->values[toMove] = 1; */
/* 	    } */

/* 	    /\* ACCEPT/REJECT STEP *\/ */
/* 	    /\* compute the likelihood *\/ */
/* 	    logRatio = loglikelihood_all(dat, dnaInfo, gen, tempPar, rng) - loglikelihood_all(dat, dnaInfo, gen, currentPar, rng); */

/* 	    /\* /\\* MH correction *\\/ *\/ */
/* 	    /\* /\\* like ratio x ( Pmove(current)/Pmove(temp) ) *\\/ *\/ */

/* 	    /\* curMut = mat_int_ij(dnaInfo->mutation1, toMove, tempPar->alpha->values[toMove]) + mat_int_ij(dnaInfo->mutation2, toMove, tempPar->alpha->values[toMove]); *\/ */
/* 	    /\* tempMut = mat_int_ij(dnaInfo->mutation1, toMove, currentPar->alpha->values[toMove]) + mat_int_ij(dnaInfo->mutation2, toMove, currentPar->alpha->values[toMove]); *\/ */

/* 	    /\* logRatio += log(mcmcPar->Pmove_alpha_old/mcmcPar->Pmove_alpha_new); *\/ */

/* 	    /\* /\\* debugging *\\/ *\/ */
/* 	    /\* ll1=loglikelihood_all(dat, dnaInfo, gen, currentPar, rng); *\/ */
/* 	    /\* ll2=loglikelihood_all(dat, dnaInfo, gen, tempPar, rng); *\/ */

/* 	    /\* if p(new/old) > 1, accept new *\/ */
/* 	    if(logRatio>=0.0) { */
/* 		/\* /\\* debugging *\\/ *\/ */
/* 		/\* printf("\naccepting automatically move from %d->%d to %d->%d (respective loglike:%f and %f)\n",vec_int_i(currentPar->alpha,toMove), toMove+1, vec_int_i(tempPar->alpha,toMove), toMove+1, ll1, ll2); *\/ */
/* 		/\* fflush(stdout); *\/ */

/* 		currentPar->alpha->values[toMove] = vec_int_i(tempPar->alpha,toMove); */
/* 		currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove); */
/* 		/\* mcmcPar->Pmove_alpha_old = mcmcPar->Pmove_alpha_new; *\/ */
/* 		mcmcPar->n_accept_alpha += 1; */
/* 	    } else { /\* else accept new with proba (new/old) *\/ */
/* 		if(log(gsl_rng_uniform(rng)) <= logRatio){ /\* accept *\/ */
/* 		    /\* /\\* debugging *\\/ *\/ */
/* 		    /\* printf("\naccepting move from %d->%d to %d->%d (respective loglike:%f and %f)\n",vec_int_i(currentPar->alpha,toMove), toMove+1, vec_int_i(tempPar->alpha,toMove), toMove+1, ll1, ll2); *\/ */
/* 		    /\* fflush(stdout); *\/ */

/* 		    currentPar->alpha->values[toMove] = vec_int_i(tempPar->alpha,toMove); */
/* 		    currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove); */
/* 		    /\* mcmcPar->Pmove_alpha_old = mcmcPar->Pmove_alpha_new; *\/ */
/* 		    mcmcPar->n_accept_alpha += 1; */
/* 		} else { /\* reject *\/ */
/* 		    tempPar->alpha->values[toMove] = vec_int_i(currentPar->alpha,toMove); */
/* 		    tempPar->kappa->values[toMove] = vec_int_i(currentPar->kappa,toMove); */
/* 		    mcmcPar->n_reject_alpha += 1; */
/* 		} */
/* 	    } /\* end  ACCEPT/REJECT STEP *\/ */
/* 	} /\* end if ancestor has changed *\/ */
/*     } /\* end for loop (for all 'i' to move) *\/ */
/* } /\* end move_alpha *\/ */






/* /\* MOVE VALUES OF KAPPA *\/ */
/* void move_kappa(param *currentPar, param *tempPar, data *dat, dna_dist *dnaInfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng){ */
/*     int i, toMove=0, temp; */
/*     double logRatio = 0.0; */


/*     /\* DETERMINE WHICH kappa_i TO MOVE *\/ */
/*     /\* sample_vec_int(mcmcPar->all_idx, mcmcPar->idx_move_kappa, FALSE, rng); *\/ */
/*     draw_vec_int_multinom(mcmcPar->all_idx, mcmcPar->idx_move_kappa, mcmcPar->move_kappa, rng); */

/*     /\* MOVE EACH kappa_i TO MOVE *\/ */
/*     for(i=0;i<mcmcPar->idx_move_kappa->length;i++){ */
/* 	toMove = vec_int_i(mcmcPar->idx_move_kappa,i); */

/* 	/\* IMPORTED CASE *\/ */
/* 	if(vec_int_i(currentPar->alpha, toMove) < 0){ */
/* 	/\* set kappa_i to 1 if imported case *\/ */
/* 	    currentPar->kappa->values[toMove] = 1; */
/* 	} else { */
/* 	    /\* INTERNAL CASE *\/ */
/* 	    /\* GET PROPOSED KAPPA_I *\/ */
/* 	    /\* needs to be on [1;maxK]*\/ */
/* 	    temp = tempPar->kappa->values[toMove] + (gsl_rng_uniform(rng) >= 0.5 ? 1 : -1); */
/* 	    if(temp < 1) { */
/* 		temp = 1; */
/* 	    } else if(temp>gen->maxK){ */
/* 		temp = gen->maxK; */
/* 	    } */
/* 	    tempPar->kappa->values[toMove] = temp; */

/* 	    /\* PROCEED TO ACCEPT/REJECT STEP ONLY IF KAPPA HAS CHANGED *\/ */
/* 	    if(vec_int_i(tempPar->kappa,toMove) != vec_int_i(currentPar->kappa,toMove)){ */
/* 		/\* MOVE ALPHA *\/ */
/* 		tempPar->alpha->values[toMove] = choose_alpha_i(toMove, dat, currentPar, mcmcPar, rng); */

/* 		/\* ACCEPT/REJECT STEP *\/ */
/* 		/\* compute the likelihood *\/ */
/* 		logRatio = loglikelihood_all(dat, dnaInfo, gen, tempPar, rng) - loglikelihood_all(dat, dnaInfo, gen, currentPar, rng); */

/* 		/\* /\\* compute the priors *\\/ *\/ */
/* 		/\* logRatio += logprior_kappa_i(toMove,tempPar) - logprior_kappa_i(toMove,currentPar); *\/ */

/* 		/\* if p(new/old) > 1, accept new *\/ */
/* 		if(logRatio>=0.0) { */
/* 		    currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove); */
/* 		    currentPar->alpha->values[toMove] = vec_int_i(tempPar->alpha,toMove); */
/* 		    mcmcPar->n_accept_kappa += 1; */
/* 		} else { /\* else accept new with proba (new/old) *\/ */
/* 		    if(log(gsl_rng_uniform(rng)) <= logRatio){ /\* accept *\/ */
/* 			currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove); */
/* 			currentPar->alpha->values[toMove] = vec_int_i(tempPar->alpha,toMove); */
/* 			mcmcPar->n_accept_kappa += 1; */
/* 		    } else { /\* reject *\/ */
/* 			tempPar->kappa->values[toMove] = vec_int_i(currentPar->kappa,toMove); */
/* 			tempPar->alpha->values[toMove] = vec_int_i(currentPar->alpha,toMove); */
/* 			mcmcPar->n_reject_kappa += 1; */
/* 		    } */
/* 		} /\* end  ACCEPT/REJECT STEP *\/ */
/* 	    } */

/* 	} /\* end if isolate is not the oldest one *\/ */
/*     } /\* end for loop (for all 'i' to move) *\/ */
/* } /\* end move_kappa *\/ */






/* /\* MOVE VALUES OF PHI *\/ */
/* void move_phi(param *currentPar, param *tempPar, data *dat, mcmc_param *mcmcPar, gsl_rng *rng){ */
/*     int i; */
/*     double logRatio=0.0; */
/*     double QCur, QTemp; */

/*     /\* GENERATE CANDIDATE VALUE FOR PHI *\/ */
/*     /\* HERE REPLACE WITH TRUNCATED LOGNORMAL (no values >1) )*\/ */
/*     do */
/*     { */
/*     	tempPar->phi = gsl_ran_lognormal(rng,log(currentPar->phi),mcmcPar->sigma_phi); */
/*     	/\* which should be the same as: *\/ */
/* 	/\* tempPar->phi = currentPar->phi*gsl_ran_lognormal(rng,0,mcmcPar->sigma_phi); *\/ */
/*     } while(tempPar->phi>1.0); */


/*     /\* ACCEPT / REJECT *\/ */
/*    /\* likelihood *\/ */
/*     logRatio += loglike_alpha_all(tempPar) - loglike_alpha_all(currentPar); */

/*     /\* prior *\/ */
/*     logRatio += logprior_phi(tempPar) - logprior_phi(currentPar); */

/*     /\* ADD CORRECTION FOR MH truncated lognormal *\/ */
/*     QCur = gsl_cdf_gaussian_P(-log(currentPar->phi),mcmcPar->sigma_phi); */
/*     QTemp = gsl_cdf_gaussian_P(-log(tempPar->phi),mcmcPar->sigma_phi); */
/*     logRatio +=  log(tempPar->phi) - log(currentPar->phi); /\* correction for lognormal *\/ */
/*     logRatio +=   log(QCur) - log(QTemp); /\* correction for truncation (no values >1) *\/ */

/*     /\* if p(new/old) > 1, accept new *\/ */
/*     if(logRatio>=0.0) { */
/* 	currentPar->phi = tempPar->phi; */
/* 	mcmcPar->n_accept_phi += 1; */
/* 	/\* printf("\nAccepting new value\n"); *\/ */
/*     } else { /\* else accept new with proba (new/old) *\/ */
/* 	if(log(gsl_rng_uniform(rng)) <= logRatio){ /\* accept *\/ */
/* 	    currentPar->phi = tempPar->phi; */
/* 	    mcmcPar->n_accept_phi += 1; */
/* 	    /\* printf("\nAccepting new value\n"); *\/ */
/* 	} else { /\* reject *\/ */
/* 	    tempPar->phi = currentPar->phi; */
/* 	    mcmcPar->n_reject_phi += 1; */
/* 	    /\* printf("\nRejecting new value\n"); *\/ */
/* 	} */
/*     } */

/* } /\* end move_phi *\/ */










/*
>>>> TESTING <<<<
*/

/* int main(){ */
/*   /\* DECLARATIONS *\/ */
/*     int TIMESPAN, i; */
/*     data *dat; */
/*     gentime *gen; */
/*     param *par, *tempPar; */
/*     dna_dist * dnaInfo; */
/*     mcmc_param * mcmcPar; */

/*     double logPrior, logLike, logPost; */

/*     /\* INITIALIZE RNG *\/ */
/*     gsl_rng *rng = create_gsl_rng(time(NULL)); */


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
/*     printf("sizes of rows in gen: "); */
/*     for(i=0;i<gen->dens->n;i++) printf("%d ", gen->dens->rows[i]->length); */


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
/*     printf("\nParameters (par)\n"); */
/*     print_param(par); */

/*     tempPar = alloc_param(3); */
/*     copy_param(par,tempPar); */
/*     printf("\nParameters (tempPar)\n"); */
/*     print_param(par); */


/*     /\* ALLOCATE MCMCPAR *\/ */
/*     mcmcPar = alloc_mcmc_param(3); */
/*     init_mcmc_param(mcmcPar, dat); */
/*     printf("\nMCMC parameters (mcmcPar)\n"); */
/*     print_mcmc_param(mcmcPar); */


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


/*     /\* MOVE MU1 *\/ */
/*     mcmcPar->sigma_mu1 = 0.001; */
/*     for(i=0;i<500;i++){ */
/* 	move_mu1(par, tempPar, dat, dnaInfo, mcmcPar, rng); */
/* 	printf("\nmu1: %.10f (reject: %d  accept: %d  ratio: %.3f)", par->mu1, mcmcPar->n_reject_mu1, mcmcPar->n_accept_mu1, (double) mcmcPar->n_reject_mu1 / mcmcPar->n_accept_mu1); */
/*     } */
/*     printf("\n"); */
/*     fflush(stdout); */

/*     /\* MOVE GAMMA *\/ */
/*     for(i=0;i<500;i++){ */
/* 	move_gamma(par, tempPar, dat, dnaInfo, mcmcPar, rng); */
/* 	printf("\ngamma: %.10f (reject: %d  accept: %d  ratio: %.3f)", par->gamma, mcmcPar->n_reject_gamma, mcmcPar->n_accept_gamma, (double) mcmcPar->n_reject_gamma / mcmcPar->n_accept_gamma); */
/*     } */
/*     printf("\n"); */
/*     fflush(stdout); */

/*     /\* MOVE TINF *\/ */
/*     for(i=0;i<500;i++){ */
/* 	move_Tinf(par, tempPar, dat, dnaInfo, gen, mcmcPar, rng); */
/* 	printf("\nTinf:"); */
/* 	print_vec_int(par->Tinf); */
/* 	printf(" (reject: %d  accept: %d  ratio: %.3f)", mcmcPar->n_reject_Tinf, mcmcPar->n_accept_Tinf, (double) mcmcPar->n_reject_Tinf / mcmcPar->n_accept_Tinf); */
/*     } */
/*     printf("\n"); */
/*     fflush(stdout); */



/*     /\* MOVE ALPHA AND KAPPA *\/ */
/*     for(i=0;i<500;i++){ */
/* 	move_alpha_kappa(par, tempPar, dat, dnaInfo, gen, mcmcPar, rng); */
/* 	printf("\nAlpha:"); */
/* 	print_vec_int(par->alpha); */
/* 	printf("\nKappa:"); */
/* 	print_vec_int(par->kappa); */
/* 	printf(" (reject: %d  accept: %d  ratio: %.3f)", mcmcPar->n_reject_alpha, mcmcPar->n_accept_alpha, (double) mcmcPar->n_reject_alpha / mcmcPar->n_accept_alpha); */
/*     } */
/*     printf("\n"); */
/*     fflush(stdout); */


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
/*     gsl_rng_free(rng); */
/*     free_data(dat); */
/*     free_gentime(gen); */
/*     free_dna_dist(dnaInfo); */
/*     free_param(par); */
/*     free_param(tempPar); */
/*     free_mcmc_param(mcmcPar); */

/*     return 0; */
/* } */





/* /\* FIND MOST LIKELY KAPPA_I (given time to infection 'T') *\/ */
/* int find_maxLike_kappa_i(int T, gentime *gen){ */
/*     int i, out=1; */
/*     double temp=0.0, currentMax=0.0; */

/*     for(i=1;i<gen->maxK;i++){ */
/* 	temp = gentime_dens(gen, T, i); */
/* 	if(currentMax < temp) { */
/* 	    currentMax = temp; */
/* 	    out = i+1; */
/* 	} */
/*     } */

/*     return out; */
/* } /\* end find_maxLike_kappa_i *\/ */



/*
  gcc instructions

  gcc -o moves matvec.c genclasses.c structures.c init.c distances.c prior.c likelihood.c moves.c -lgsl -lgslcblas -Wall -g


  gcc -o moves matvec.c genclasses.c structures.c init.c distances.c prior.c likelihood.c moves.c -lgsl -lgslcblas -Wall -O3

 ./moves

  valgrind --leak-check=full -v moves

*/


/***************************************************************************
                            parameters.h
                         -------------------
    last modified   : 1/5/2008
    copyright       : (C) 2008 by Dhananjay Thiruvady
    libraries		: .
    description		: paramters to be used with bacs_smjs
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef PARAMETERS_H
#define PARAMETERS_H

//other definitions
#define MAX_VAL 100000			/* Global best cost*/
#define var_size 31

//learning rates
#define alp 1.0
#define bet 6.0
//#define q_0 0.9
#define rho 0.01
//#define lrate 0.01

// tau_min, tau_max: lower, respecively upper, limit for the pheromone values
#define tau_min 0.001
#define tau_max 0.999

// aco versions
#define acs 1 
#define maxmin 2
//#define aco_version 2

// solution type, plain aco, beam_aco, cpaco
#define aco 1
#define beam_aco 2
#define cp_aco 3

// heuristics
#define edd 1
#define sst 2
//#define type 3 // this is what the code will use

// dependent model
#define dependent 0

// samples for stochastic sampling
#define nsamples 5

// maximum children to generate
#define max_children 20

// how to select from children
#define heu_select 1
#define uni_select 2

// samples for beta sampling
#define n_beta_samples 100

// which of the above to use? -> heuristic
//#define selection 1

#endif

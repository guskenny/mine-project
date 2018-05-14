/******************************************************************************** 
                            Solver 
                         -------------------------- 
    last modified   : 13/5/2016 
    copyright       : (C) 2016 by Dhananjay Thiruvady 
    libraries		: . 
    description		: contains a data structure for the Solver class
********************************************************************************/ 
 
/*************************************************************************** 
 *                                                                         * 
 *   This program is free software; you can redistribute it and/or modify  * 
 *   it under the terms of the GNU General Public License as published by  * 
 *   the Free Software Foundation; either version 2 of the License, or     * 
 *   (at your option) any later version.                                   * 
 *                                                                         * 
 ***************************************************************************/ 

#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <new>
#include <vector>
#include <utility>
#include <list>

using namespace std;

#include "../include/daten.h"
#include "MaxClosure_NetworkFlow_LR.h"
#include "Timer.h"
#include "mip.h"
#include "solver_functions.h"
#include "../include/graph.h"
#include "../include/BranchNode.h"
#include "ACO_solution.h"
#include "pheromones.h"
#include "../parallel/MineProblem.h"

/*
  A helper class: Solver functions
  Contains functions that solver needs
*/

/*
  The Solver class 
  Esentially, every algorithm based on a network 
    will sit within this class
 */

class Solver{

 private:
  Daten *data;
  SolverFunctions *funcs;

  // some local functions
  vector<vector<double> > getProfit(const vector<vector<double> > &lambda);
  vector<vector<int> > createPreds();
  double Average(const vector<double> &run_times);
  // algorithm related functions
  void Primal(MaxClosure_NetworkFlow_LR *net, const vector<double> &profit);
  void Dual(const vector<vector<double> > &lambda, 
	    double &upper_bound, bool show);   // needs to be implemented
  void buildNetwork(MaxClosure_NetworkFlow_LR *net);
  vector<vector<vector<double> > > runNetworkFlow(double &upper_bound, bool show, MaxClosure_NetworkFlow_LR *net);
  void MIP_run(bool simple);
  void NF_run(int type);  
  ACO_Solution* ACO(Pheromones *p, int nants, Random *r, double q_0, double lrate);

 public:
  Solver(Daten &data, SolverFunctions &funcs);
  // Run an ACO with heuristics
  ACO_Solution* ACO_run(Pheromones *p, ACO_Solution *lr_sol, int max_time, int nants, double lrate, double q_0, bool disp, bool t_limit);
  ACO_Solution* ACO_run(Pheromones *p, ACO_Solution * lr_sol, int max_time, 
			int nants, double lrate, double q_0, bool disp, 
			bool t_limit, MineProblem *problem);
  // Run the large neighbourhood search
  double LargeNeighbourhoodSearch(bool disp, 
				  int time_limit, bool show, ACO_Solution *sol,
				  vector<int> &X, 
				  vector<vector<double> > &Y,
				  MineProblem *problem 
				  );

  // Run a window search
  double WindowSearch( bool disp, 
		       int time_limit, bool show, ACO_Solution *sol,
		       vector<int> &X, 
		       vector<vector<double> > &Y,
		       MineProblem *problem,
		       int time, int shift
		       );
    
  // Like the above, but fixed when the blocks can't be done via BranchNode_info
  double LNS_fixed(bool disp, 
		   int time_limit, bool show,
		   MineProblem *problem,
		   vector<int> &X, 
		   vector<vector<double> > &Y 
		   );
  
  // Run a Lagrangian relaxation
  ACO_Solution* Lagrangian_heuristic(bool disp, double time_limit, bool show, 
				     const double &bound, int iterations, int nants, double lrate, 
				     double q_0, bool use_aco, MineProblem *problem);
    
  // Or specify which algorithm to run through a shell
  bool run_algorithm(bool disp, int time_limit, bool show, int type, 
		     bool simple, int nf_type, double bound, int nants, double lrate, 
		     double q_0, bool use_aco);



  ~Solver();
};


#endif

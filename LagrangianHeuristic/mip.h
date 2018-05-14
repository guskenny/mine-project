/*************************************************************************** 
                            MIP algorithms 
                         ------------------- 
    last modified   : 13/5/2016 
    copyright       : (C) 2016 by Dhananjay Thiruvady 
    libraries		: . 
    description		: contains all the MIP algorithms and 
                          required functions 
 ***************************************************************************/ 
 
/*************************************************************************** 
 *                                                                         * 
 *   This program is free software; you can redistribute it and/or modify  * 
 *   it under the terms of the GNU General Public License as published by  * 
 *   the Free Software Foundation; either version 2 of the License, or     * 
 *   (at your option) any later version.                                   * 
 *                                                                         * 
 ***************************************************************************/ 

#ifndef MIP_H
#define MIP_H

#include <iostream>
#include <new>
#include <vector>
#include <string>
#include <array>
#include "gurobi_c++.h"

using namespace std;

#include "Random.h"
#include "../include/daten.h"
#include "../include/BranchNode.h"
#include "Timer.h"
#include "../parallel/MineProblem.h"
#include "ACO_solution.h"
#include "solver_functions.h"

class MIP{

 private:
  Daten *data;

  // some local functions

  // Create the variables in the model
  void setVariables( GRBModel &model, 
		     vector<vector<GRBVar> > &Xvar,  
		     vector<vector<vector<GRBVar> > > &Yvar);
    
  // Set the objective
  void setObjective( GRBModel &model, 
		     vector<vector<vector<GRBVar> > > &Yvar);

  // Precedence constraints
  void precedenceConstraints( GRBModel &model, 
			      vector<vector<GRBVar> > &Xvar);

  // Precedence constraints for the extended model
  void precedenceConstraintsExt( GRBModel &model, 
				 vector<vector<GRBVar> > &Xvar);


  // Sum destination and Clique block constraints
  void sumDestCliqueBlockConstraint(GRBModel &model,
				     vector<vector<GRBVar> > &Xvar,
				     vector<vector<vector<GRBVar> > > &Yvar);

  
  //  Sum destination and block completion constraints
  //  For the extended model
  void sumDestCompletionConstraint( GRBModel &model, 
				    vector<vector<GRBVar> > &Xvar,
				    vector<vector<vector<GRBVar> > > &Yvar);
    
  // Resource constraints
  void resourceConstraints( GRBModel &model, 
			    vector<vector<vector<GRBVar> > > &Yvar);
  
  // Other functions
  void resetY(vector<vector<vector<GRBVar> > > &Yvar);
  bool isResourceFeasible(const vector<vector<vector<double> > > &Yvar_sol);
  double getNegAvgProfit(const vector<double> &profit);
  void updateSolution(const vector<vector<vector<double> > > &Yvar_sol, 
		      vector<int> &X, vector<vector<double> > &Y);
  void updateParSol(MineProblem *problem, const double &obj,
		    vector<vector<vector<double> > > &Yvar_sol
		    );
  void generateParSolution(const vector<vector<vector<double> > > &Yvar_sol,
			   const double &tobj, Sol_Int &sol);
  void setModelParams(GRBModel &model);
  void checkBlockAssignments(const vector<vector<vector<double> > > &solution, 
			     const vector<vector<vector<double> > > &Yvar_sol,
			     const vector<double> &profit);
  void updateMaster(MineProblem *problem, const double &obj,
		    const vector<vector<vector<double> > > &Yvar_sol);

  double checkACOSolution( const vector<vector<vector<double> > > &Yvar_sol);

 public:
   MIP(Daten &data);
   void MIP_simple();
   void MIP_extended();
   vector<vector<vector<double> > > MIP_relaxed(const vector<vector<double> > &lambda, 
						double &upper_bound, bool show);
   vector<vector<vector<double> > > MIP_relaxed_efficient(const vector<vector<double> > &lambda, 
							  double &upper_bound, bool show);

   void MIP_fixed(GRBModel &model, vector<vector<vector<GRBVar> > > &Yvar);
   // Set which blocks can be mined anywhere
   // and solve
   double MIP_fixed_resolve(int mod_tmin, int mod_tmax, 
			    bool show, GRBModel &model, vector<vector<vector<GRBVar> > > &Yvar,
			    vector<vector<vector<double> > > &Yvar_sol,
			    const vector<double> &profit
			    );

   // For the variant where you only allow a period where the block can be mined
   double MIP_fixed_resolve(bool show, GRBModel &model, vector<vector<vector<GRBVar> > > &Yvar,
			    vector<vector<vector<double> > > &Yvar_sol,
			    MineProblem *problem
			    );

   double VLNS(const vector<vector<vector<double> > > &solution, 
	       const vector<double> &profit, Timer &timer,
	       vector<int> &X,
	       vector<vector<double> > &Y,
	       MineProblem *problem
	       );
   double VLNS_BW(const vector<vector<vector<double> > > &solution, 
		  const vector<double> &profit, Timer &timer,
		  vector<int> &X,
		  vector<vector<double> > &Y,
		  MineProblem *problem
		  );
   double MIP_fixed_partial(const vector<bool> &blku, int mod_tmin, int mod_tmax, 
			    vector<vector<vector<double> > > &Yvar_sol,
			    const vector<vector<double> > &X,
			    bool last_window
			    );

   double VLNS_Small(const vector<vector<vector<double> > > &solution, 
		     const vector<double> &profit, Timer &timer,
		     vector<int> &X,
		     vector<vector<double> > &Y,
		     MineProblem *problem
		     );
   double VLNS_fixed(const vector<double> &profit, Timer &timer,
		     vector<int> &X,
		     vector<vector<double> > &Y,
		     MineProblem *problem
		     );

   double VLNS_Window(const vector<vector<vector<double> > > &solution, 
		      const vector<double> &profit, Timer &timer,
		      vector<int> &X,
		      vector<vector<double> > &Y,
		      MineProblem *problem,
		      int time, int shift
		      );
   ~MIP();
};

#endif

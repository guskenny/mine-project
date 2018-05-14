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

#ifndef SOLVER_FUNCS_H
#define SOLVER_FUNCS_H

#include <iostream>
#include <new>
#include <vector>
#include <set>
#include <utility>

using namespace std;

#include "../include/daten.h"
#include "MaxClosure_NetworkFlow_LR.h"
#include "../include/BranchNode.h"
#include "ACO_solution.h"

/*
  A helper class: Solver functions
  Contains functions that solver needs
*/


class SolverFunctions{

 private:
  Daten *data;
  vector<vector<int> > max_block_profit;
  int diff_amt;

  // some local functions
  double ComputeAvgResourceUse(int r);
  // 
  void setProfit();
  void setBackPropProfit();
  void setBackPropProfitAvg();
  void setPosBackPropProfit();
  void setSucc();
  void setNegBlock();
  bool getNegSucc(int block, vector<pair<bool, bool> > &mem_succ);
  double getNegSucc(int block, vector<pair<double, bool> > &mem_succ);
  double getSuccProfit(int block, vector<pair<double, bool> > &mem_profit);
  void dispSucc(int block, vector<bool> &visited);
 public:
  // data for ACO solutions
  vector<double> profit;
  vector<int> block_dest;
  vector<set<int> > succ;
  vector<bool> negBlocks;

  SolverFunctions(Daten &data, int diff_amt);
  vector<vector<double> > update_lambda(double gap, const vector<vector<double> > &lambda, 
					double delta, const vector<vector<vector<double> > > &Yvar);
  double ComputeAvgBlockProfitSum();  
  vector<vector<vector<double> > > generateSolution(MaxClosure_NetworkFlow_LR *net);

  // profit related functions
  void setProfit(  vector<double> &profit,
		   vector<vector<double> > &profit_time);
  void setLambdaProfit(  vector<vector<double> > &profit_time,
			 const vector<vector<double> > &lambda);
  void setMaxBlockProfit();
  void updateProfit(const vector<vector<double> > &lambda);
  void setLambda(vector<vector<double> > &lambda);
  bool validateSols(const vector<vector<vector<double> > > &Yvar,
		    const vector<vector<vector<double> > > &Yvar_mod);
  void setDeepBackPropProfit();

  void generateParSolution(ACO_Solution *a, Sol_Int &sol);

  ~SolverFunctions(){};
};


#endif

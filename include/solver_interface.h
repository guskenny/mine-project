

/******************************************************************************** 
                            Solver Interface
                         -------------------------- 
    last modified   : 25/08/2016 
    copyright       : (C) 2016 by Dhananjay Thiruvady 
    libraries: . 
    description: contains an iterface to the solvers
********************************************************************************/ 
 
/*************************************************************************** 
 *                                                                         * 
 *   This program is free software; you can redistribute it and/or modify  * 
 *   it under the terms of the GNU General Public License as published by  * 
 *   the Free Software Foundation; either version 2 of the License, or     * 
 *   (at your option) any later version.                                   * 
 *                                                                         * 
 ***************************************************************************/ 

#ifndef SOLVER_INT_H
#define SOLVER_INT_H

#include <iostream>
#include <vector>

using namespace std;

#include "../include/CumulativeModel.h"
#include "../LagrangianHeuristic/solver_functions.h"
#include "../LagrangianHeuristic/solver.h"
#include "../LagrangianHeuristic/ACO_solution.h"
#include "../parallel/MineProblem.h"
#include "../include/BranchNode.h"
/*
  The Solver Interface class 
  Allows creating objects of different solvers
  
*/

class SolverInterface{

 protected:
  CumulativeModel *data;
  bool deleteData;
  Sol_Int * mip_sol;

  // Some local functions

  // Generate a Y variable like input for the algorithms
  void convertToY(const BranchNode_info &bni, vector<vector<vector<double> > > &Y);
  // Convert the current solution to an integer solution
  Sol_Int* convertToSol_Int(const vector<int> &X, 
			    const vector<vector<double> > &Y, 
			    const double &obj
			    );
  void convertToSol_Int(vector<int> &X, vector<vector<double> > &Y);
  // collect data from ACO to store in problem
  void collecDataFromACO(ACO_Solution *a, vector<int> &X, vector<vector<double> > &Y);
  // collect data from single cumulative real vector and store in sol_real
  // x/y vectors below must be length data.graph.getNumNodes()
  void storeCumulativeSoln(Sol_Real &sol,const vector<double> &y);
  void storeCumulativeSoln(Sol_Int &sol,const vector<double> &y); 
  void storeCumulativeSoln(Sol_Int &sol,const vector<int> &x);

 public:
  SolverInterface(Daten &data);
  SolverInterface(CumulativeModel &data) { this->data = &data; deleteData= false;}
  SolverInterface(SolverInterface &other) ;
  
  // Call one of the solvers (e.g. Lagrangian relaxation) using the Y variable
  // Returns success (1) or failure (0), or maybe something else?
  // possible alg_type values:
  // 1   ACO
  // 2   BZ
  // 3   LR+ACO
  // 4 VLNS
  // 5   VLNS with   2 time periods
  // 6   Volume Algorithm
  // 7   LaPSO  
  // char(alg_type) == 'L' LaPSO
  int callSolver(MineProblem &problem, int diff_amt, int iterations, int alg_type, int time, int shift);
    
  // call solver with default arguments/options   
  virtual  int callSolver(MineProblem &problem) {
    return callSolver(problem,0,1,2,-1,-1); // use BZ
  }

  // Update BNI with the latest solution information
  void updateBNI(MineProblem &problem);

  // Update a solution, required by the main parallel solver
  void updateSolution(Sol_Int &sol); // THIS FUNCTION DOES NOTHING!!!!

  // Return a solution determined by the MIP
  Sol_Int *getMipSol(){
    return mip_sol;
  } 
  // Maybe delete data
  virtual ~SolverInterface();
};


#endif


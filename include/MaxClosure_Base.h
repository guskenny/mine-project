/*************************************************************************** 
                            MaxClosure_Base 
                         ------------------- 
    last modified   : 13/5/2016 
    copyright       : (C) 2016 by Dhananjay Thiruvady 
    libraries		: . 
    description		: contains a data structure for the MaxClosure class
 ***************************************************************************/ 
 
/*************************************************************************** 
 *                                                                         * 
 *   This program is free software; you can redistribute it and/or modify  * 
 *   it under the terms of the GNU General Public License as published by  * 
 *   the Free Software Foundation; either version 2 of the License, or     * 
 *   (at your option) any later version.                                   * 
 *                                                                         * 
 ***************************************************************************/ 

#ifndef MCBase_H
#define MCBase_H

#include <iostream>
#include <new>
#include <vector>

#include "graph.h"

/*
  The Max Closure class
 */

class MaxClosure_Base{

protected:
  int _status;
  const Graph *_graph;

 public:
  //MaxClosure_Base(const char * filename);
  MaxClosure_Base(const Graph &graph);
  virtual ~MaxClosure_Base();
  
  int getStatus();

  virtual void setProfit(const std::vector<double> &profit)=0; // profit for each vertex of the network
  virtual int solve() = 0; 	// returns non-zero if an error occured
  virtual void getResidualProfit(std::vector<double> &residualProfit)=0;
  // getClosure returns true if a solution was found
  virtual int getClosure(std::vector<int>  &sol) = 0; // 0/1 for each vertex;

  virtual double calcProfit(const std::vector<double> &profit,const std::vector<int> &soln) const;
  // validation of a solution: return false if solution doesn't satisfy precedences
  virtual bool isClosure(const std::vector<int> &soln) const; 


  // get dual solution information:
  //virtual bool getInFlowSolution(std::vector<std::vector<double> > &inFlow)=0;  //inFlow[i][j] is flow from node pred[i][j] to i
  // outFlow is really the same information as inFlow
  // virtual bool getOutFlowSolution(std::vector< std::vector<double> > &outFlow)=0; //outFlow[i][j] is flow from node i to succ[i][j] to i
  //virtual bool getRedCost(std::vector<double> &sol) = 0; // reduced cost (dual) associated with 0<= x_i <= 1 constraints
  //virtual void MIP_run(const char *filename, bool simple)=0;
  //virtual void NF_run(const char *filename, int type)=0;  
  //virtual void Lagrangian_heuristic(const char *filename, bool disp, double time_limit, bool show)=0;
  

};
#endif

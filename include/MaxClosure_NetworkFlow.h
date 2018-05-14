/*************************************************************************** 
                            A Network Class 
                         ------------------- 
    last modified   : 13/5/2016 
    copyright       : (C) 2016 by Dhananjay Thiruvady 
    libraries		: . 
    description		: contains a data structure for a network
 ***************************************************************************/ 
 
/*************************************************************************** 
 *                                                                         * 
 *   This program is free software; you can redistribute it and/or modify  * 
 *   it under the terms of the GNU General Public License as published by  * 
 *   the Free Software Foundation; either version 2 of the License, or     * 
 *   (at your option) any later version.                                   * 
 *                                                                         * 
 ***************************************************************************/ 

#ifndef MCNF_H
#define MCNF_H

#include <ilcplex/ilocplex.h> 
#include <iostream>
#include <new>
#include <vector>
#include <string>
#include <set>
#include "graph.h"
#include "MaxClosure_Base.h"

/*
  A structure to pass around CPXNET data
 */

struct CPXData{
  int nnodes;
  int nedges;
  // associated with edges
  std::vector<int> from;
  std::vector<int> to; 
  std::vector<double> low;
  std::vector<double> up; 
  std::vector<double> obj;
  // associated with nodes
  std::vector<double> supply;
  void reserve(int nn,int ne){
		nnodes=nedges=0;
		from.reserve(ne);
		to.reserve(ne);
		low.reserve(ne);
		up.reserve(ne);
		obj.reserve(ne);
		supply.reserve(nn);
  }
};

/*
  Minimalist implementation of MaxClosure solution via network flow
*/

class MaxClosure_NetworkFlow : public MaxClosure_Base {

 private:
  // dummy node, may or may not be used
//  const int dummyNd=0;
//  const int offset=1; // allows for dummy node
//  const double scale=1e6; // scale factor for supply/flows
  const int dummyNd;
  const int offset; // allows for dummy node
  const double scale; // scale factor for supply/flows
  // some CPLEX Network data
  CPXENVptr env; // may need to be static for parallel running
  CPXNETptr prob;

 public:
  MaxClosure_NetworkFlow(const Graph &graph);
  ~MaxClosure_NetworkFlow();
	  
  // cplex network related functions
  int updateSupply();

  // inherited from MaxClosureBase
  virtual int solve();
  virtual void setProfit(const std::vector<double> &profit);
  virtual int getClosure(std::vector<int>  &sol);
  void getResidualProfit(std::vector<double> &residualProfit){}
  virtual double getObjective() const;
};


#endif

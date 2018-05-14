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

#ifndef MCNF_LR_H
#define MCNF_LR_H

#include <ilcplex/ilocplex.h> 
#include <iostream>
#include <new>
#include <vector>
#include <string>
#include <set>

using namespace std;

#include "../include/graph.h"
#include "../include/lr_graph.h"
#include "../include/MaxClosure_Base.h"
#include "../include/BranchNode.h"


/*
  A structure to pass around CPXNET data
 */

struct CPXData_LR{
  int nnodes;
  int nedges;
  // associated with edges
  vector<int> from;
  vector<int> to; 
  vector<double> low;
  vector<double> up; 
  vector<double> obj;
  // associated with nodes
  vector<double> supply;
  vector<int> nodeIds;
  vector<int> edgeIds;
};

/*
  A generic Network class 
  Supposed to represent network that is usable with 
    any method, LR or MC, etc.
 */

class MaxClosure_NetworkFlow_LR : public MaxClosure_Base {

 private:
  // A dummy node, may or may not be used
  NetExtNode *dummy;
  // Block by Time representation
  vector<vector<NetExtNode*> > nodeVarMatrix;
  // An expansion of the above into a single vector for 
  //   the purposes CPLEX's NetowkrFlow Algorithm
  // Dummy node to be included here
  vector<NetExtNode*> nodeVar; 
  // Edges
  vector<NetExtArc*> edgeVar;
  // mapping for nodes
  vector<pair<int,int> > nodeIdMapping;
  // variables related to profits, etx
  double profitSum;
  vector<vector<double> > profit; // always associated with time
  // some CPLEX Network data
  CPXENVptr env;
  CPXNETptr prob;

 public:
  //MaxClosure_NetworkFlow_LR(const vector<vector<double> > &profit, Graph &graph);
  MaxClosure_NetworkFlow_LR(const vector<vector<double> > &profit, Graph &graph, 
			    const vector<vector<int> > &preds,
			    const BranchNode_info *bni
			    ); 
  
  // cplex network related functions
  int updateSupply();
  int solve();
  void setProfit(const vector<double> &profit);
  int getClosure(vector<int>  &sol);

  double getSolInfo();
  // network related functions
  void setExtNodeVal(int nodeId, double value);
  double getExtNodeVal(int block, int time);
  void getSupplyInfo(pair<vector<int>, vector<double> > &supply);
  CPXData_LR* getExtNodeExtArcInfo();
  void resetValues();
  void updateProfit(const vector<vector<double> > &profit);
  // some index related functions
  int getProfitSize(){ return profit.size();};
  void setProfitSum();
  ~MaxClosure_NetworkFlow_LR();

};


#endif

/*************************************************************************** 
                     A template class for boost max flow
                    -------------------------------------
    last modified   : 21/6/2016
    copyright       : (C) 2016 by Angus Kenny
    libraries		    : .
    description		  : contains data structure for boost algorithms, inherits
                      from MaxClosure_Base class. Class is a template class
                      which allows different boost algorithms to be used.

                      Use the classes:

                      MaxClosure_BoostMaxFlow_PL(graph) for push-relabel
                      MaxClosure_BoostMaxFlow_BK(graph) for Boykov-Kolmogorov
                      MaxClosure_BoostMaxFlow_EK(graph) for Edmonds-Karp

                      Function definitions have not been separated from the
                      header file because it is a template.

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef PreprocMC_H
#define PreprocMC_H

#include "MaxClosure_Base.h"
#include <map>
#include <set>
#include "MaxClosureFactory.h"
#include "boostMaxFlow.h"

struct PP_Node {
    //Graph * graph; 			 // pointer to reduced graph
  int nodeID; 	 			 // value in reduced graph
  std::vector<int> * chain;  // pointer to chain (if in chain)
  int status;  				 // status (either 0, 1 or -1 if not fixed and -2 if chain)
    PP_Node(Graph &_graph) {}//: graph(&_graph){}
};

class MaxClosure_PP : public MaxClosure_Base {
  private:
    int reducedNumNodes;
    std::vector< std::vector<int> > chains; // vector of chains
    std::vector<PP_Node> pp_nodes;
    MaxClosure_Base * reducedProblem;
    Graph reducedGraph;

  public:
    MaxClosure_PP(const Graph &graph, const std::map <int,int> &fixed, MaxClosureFactory mcfactory);
    ~MaxClosure_PP();

    // Virtual function, to be implemented by specific class
    int solve();
    void setProfit(const std::vector<double> &profit); // profit for each vertex of the network
    int getClosure(std::vector<int> &sol); // 0/1 for each vertex
    void getResidualProfit(std::vector<double> &residualProfit);
    void visitDFS_PP(int curr, std::vector<int> &status, std::vector<int> &inDegree, std::vector<int> &outDegree, int value, bool backwards);
    void findChain(NodeID curr, std::vector<int> &inDegree, std::vector<int> &outDegree, std::vector<int> &status, std::vector<int> &chain, bool backwards);

    // get dual solution information:
    //bool getInFlowSolution(std::vector< std::vector<double> > &inFlow){};
    //inFlow[i][j] is flow from node pred[i][j] to i
    // outFlow is really the same information as inFlow
    // virtual bool getOutFlowSolution(std::vector< std::vector<double> > &outFlow)=0;
    //outFlow[i][j] is flow from node i to succ[i][j] to i
    //bool getRedCost(std::vector<double> &sol){};
    // reduced cost (dual) associated with 0<= x_i <= 1 constraints
};

#endif

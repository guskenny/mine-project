/*************************************************************************** 
                          Implementation of Network class 
                         ------------------------------------------- 
    last modified   : 13/5/2016 
    copyright       : Monash University 2016
    libraries		: . 
    description		: the network data structure
 ***************************************************************************/ 
 
/*************************************************************************** 
 *                                                                         * 
 *   This program is free software; you can redistribute it and/or modify  * 
 *   it under the terms of the GNU General Public License as published by  * 
 *   the Free Software Foundation; either version 2 of the License, or     * 
 *   (at your option) any later version.                                   * 
 *                                                                         * 
 ***************************************************************************/ 

#include "MaxClosure_NetworkFlow.h" 

/*
  Set up the network in the constructor
*/

MaxClosure_NetworkFlow::MaxClosure_NetworkFlow(const Graph &graph) : MaxClosure_Base(graph), dummyNd(0), offset(1), scale(1e6) {

  // set up some parameters
  const int nn = graph.getNumNodes()+offset;
  int status;
  try{
    env = CPXopenCPLEX (&status);
    if(status){
      std::cerr << "Error opening environment ... " << status << std::endl;
      _status = status;
      return;
    }
    //CPXsetintparam(env,CPX_PARAM_SCRIND,1); // useful during debugging
    prob = CPXNETcreateprob(env, &status, "NetFlowProb");
    if(status){
      std::cerr << "Error creating network flow problem... " << status << std::endl;
      _status = status;
      return;
    }
    CPXData cd;
    cd.reserve(nn,graph.getNumArcs()+nn/2); // allow for arcs to/form dummy
    cd.supply.resize(nn,0.0);
    // set up edges based on graph:
    for(int a=0;a<graph.getNumArcs();++a){
      const Arc *arc=graph.getArc(a);
      cd.from.push_back(arc->getSrcID()+offset);
      if( cd.from.back() >= nn || cd.from.back() < 0)
	std::cerr << "From node idx out of range " << cd.from.back() << std::endl;
      cd.to.push_back(arc->getTgtID()+offset);
      if( cd.to.back() >= nn || cd.to.back() < 0)
	std::cerr << "To node idx out of range " << cd.to.back() << std::endl;
    }
    cd.obj.resize(graph.getNumArcs(),0.0);
    for(int n=0;n<graph.getNumNodes();++n){
      const Node *node=graph.getNode(n);
      if(node->getInDegree() == 0){
	cd.from.push_back(dummyNd);
	if( dummyNd >= nn || dummyNd < 0)
	  std::cerr << "Dummy from idx out of range " << dummyNd << std::endl;
	cd.to.push_back(node->getID()+offset);
	if( node->getID()+offset >= nn || node->getID()+offset < 0)
	  std::cerr << "To idx out of range " << (node->getID()+offset) << std::endl;
	cd.obj.push_back(1.0);
      }else if(node->getOutDegree() == 0){
	cd.from.push_back(node->getID()+offset);
	if( node->getID()+offset >= nn || node->getID()+offset < 0)
	  std::cerr << "From idx out of range " << (node->getID()+offset) << std::endl;
	cd.to.push_back(dummyNd);
	if( dummyNd >= nn || dummyNd < 0)
	  std::cerr << "Dummy to idx out of range " << dummyNd << std::endl;
	cd.obj.push_back(0.0);
      }
    }
    int ne=(int)cd.obj.size();
    cd.low.resize(ne,0.0);
    cd.up.resize(ne,CPX_INFBOUND); // infinite capacity
    int objsen = CPX_MIN;
    status = CPXNETcopynet(env, prob, objsen,
			   nn, &cd.supply[0], NULL, // no node names
			   ne, &cd.from[0], &cd.to[0], &cd.low[0],
			   &cd.up[0], &cd.obj[0], NULL); // no arc names
    if(status){
      std::cerr << "Error creating CPXNET ... " << status << std::endl;
      _status = status;
      return;
    }
  }catch(IloException &e){
    std::cerr << "Error in cplex network construction ..." << e<< std::endl;
  }catch(...){
    std::cerr << "Unknown exception during cplex network construction ... " << std::endl;
  }
} // end MaxClosure_NetworkFlow constructor


/*
  Clean up the memory allocated
*/
MaxClosure_NetworkFlow::~MaxClosure_NetworkFlow(){

  CPXNETfreeprob(env,&prob);
  prob=0;
  CPXcloseCPLEX(&env);
}




/*
  Optimise
    Objective must be scaled after
*/

int MaxClosure_NetworkFlow::solve(){
  _status= CPXNETprimopt(env,prob);
  int primFeas,dualFeas;
  CPXNETsolninfo(env,prob,&primFeas,&dualFeas);
//  std::cout <<"Network solution is primal "
//	    << (primFeas ? "feasible" : "infeasible")
//	    << " and dual "
//	    << (dualFeas ? "optimal" : "infeasible")
//	    <<std::endl;
  return _status;
}


/*
  Return the node values after a solve
*/

int MaxClosure_NetworkFlow::getClosure(std::vector<int>  &sol){
  int nn=CPXNETgetnumnodes(env,prob);
  std::vector<double> pi(nn);
  int status = CPXNETgetpi(env, prob, &pi[0], 0, nn-1);
  if(status){
    std::cerr << "Error obtaining network solution (node values) ... " << status << std::endl;
    _status = status;
    return false;
  }
  sol.resize(nn-1);	// make sure it's the right size
  const double eps=1e-3; // can be large as we are truncating
  for(int i=0;i<sol.size();++i)
    sol[i] = int(fabs(pi[i+offset] - pi[dummyNd])+eps);
  return true;
}

/*
  Set the profit for each node in the network
  Update supply values without re-creating the network
  Notes:
    The vector "profit" must be the same length as number of vertices (numnodes-1)
    The supplies are scaled: val/1e8, So the objective needs to be scaled appropriately
*/
void MaxClosure_NetworkFlow::setProfit(const std::vector<double> &profit){
  int nn=CPXNETgetnumnodes(env,prob);
  std::vector<double> supply(nn);
  std::vector<int> nodeIds(nn);
  for(int i=0;i<nn;++i) nodeIds[i] = i;
  supply[dummyNd]=0;
  for(int i = 0; i < (int)profit.size(); i++){
    supply[i+offset] = -profit[i]/scale;
    supply[dummyNd] -= supply[i+offset];
  }
  _status = CPXNETchgsupply (env, prob, nn, &nodeIds[0], &supply[0]);
}

double MaxClosure_NetworkFlow::getObjective() const {
  double obj;
  CPXNETgetobjval(env,prob,&obj);
  return obj*scale;
}



/*************************************************************************** 
                          Implementation of Network class 
                         ------------------------------------------- 
    last modified   : 13/5/2016 
    copyright       : (C) 2016 by Dhananjay Thiruvady 
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

#include <ilcplex/ilocplex.h> 

using namespace std; 

#include "MaxClosure_NetworkFlow_LR.h" 


ILOSTLBEGIN

/*
  Set up the network in the constructor
*/

MaxClosure_NetworkFlow_LR::MaxClosure_NetworkFlow_LR(const vector<vector<double> > &profit, Graph &graph, 
						     const vector<vector<int> > &preds,
						     const BranchNode_info *bni
						     ) : 
  MaxClosure_Base(graph){

  // set up some parameters
  const int nB = profit.size();
  const int t_max = profit[0].size();

  this->profit = profit;
  // Set up network
  // A dummy node to conserve the flow
  // The first ID is 0
  setProfitSum();
  dummy = new NetExtNode(0, profitSum, 0.0); 
  // the first node has no mapping to meaningful blocks and times
  nodeIdMapping.push_back(make_pair(-1,-1)); 

  nodeVar.reserve(nB);
  //edgeVar.reserve(nB*nB);

  // Set up all the nodes
  nodeVarMatrix.resize(nB);
  int count = 0;
  for(int b = 0; b < nB; b++ ){
    nodeVarMatrix[b].resize(t_max);
    for(int t = 0; t < t_max; t++ ){
      count++; // should start at 1 due to the dummy
      nodeVarMatrix[b][t] = new NetExtNode(count, -profit[b][t], 0.0);
      nodeIdMapping.push_back(make_pair(b,t));
    }
  }
  // set up the nodes in a single vector
  vector<int> succ_cost(nB,0);
  
  // Set up the arcs, based on precedences
  // Consider time as a precedence as well
  // First within a timeslot
  // New addition (27-9-16)
  // Include bni information when creating nodes
  for(int a=0; a<nB; a++){
    //    std::vector<int> pred = graph.getPreds(a);
    std::vector<int> pred = preds[a];
    int n = pred.size();
    for(int p=0; p<n; p++){
      int b = pred[p];
      for(int t = 0; t < t_max; t++ ){
	//if(bni.time[b][0] != bni.time[b][1] &&
	//   t >= bni.time[b][0] && t <= bni.time[b][1]-1)
	  edgeVar.push_back(new NetExtArc(edgeVar.size(), nodeVarMatrix[b][t], nodeVarMatrix[a][t], 0, 1e99, 0));
	  //else // the edges outside this time can't be used
	      // edgeVar.push_back(new NetExtArc(edgeVar.size(), nodeVarMatrix[b][t], nodeVarMatrix[a][t], 0, 0, 0));
      }
    }
  }
  // Now within a block
  for(int t = 0; t < t_max-1; t++ ){
    for(int b=0; b<nB; b++){
      //if(bni.time[b][0] != bni.time[b][1] 
      // && t >= bni.time[b][0] && t <= bni.time[b][1]-1)
	edgeVar.push_back(new NetExtArc(edgeVar.size(), nodeVarMatrix[b][t+1], nodeVarMatrix[b][t], 0, 1e99, 0));
	//else // the edges outside this time can't be used
	//edgeVar.push_back(new NetExtArc(edgeVar.size(), nodeVarMatrix[b][t+1], nodeVarMatrix[b][t], 0, 0, 0));
    }
  }
  //cout << "Predecessors: " << edgeVar.size() << endl;

  count = 0;
  // Set up the edges for the dummy node
  for(int b=0; b<nB; b++){
    //int n = graph.getPreds(b).size();
    int n = preds[b].size();
    if( n == 0 ){ // we have no predecessors
      count++;
      edgeVar.push_back(new NetExtArc(edgeVar.size(), dummy, nodeVarMatrix[b][t_max-1], 0, 1e99, 1.0));
    }
    int n1 = succ_cost[b];
    if( n1 == 0 ){ // we have no predecessors
      edgeVar.push_back(new NetExtArc(edgeVar.size(), nodeVarMatrix[b][0], dummy, 0, 1e99, 0));
    }
  }
  //cout << "No precs: " << count  << " and total edges: " << edgeVar.size() << endl;

  // spread out the matrix into nodeVar
  nodeVar.push_back(dummy);
  for(int b = 0; b < nB; b++ ){
    for(int t = 0; t < t_max; t++ ){
      nodeVar.push_back(nodeVarMatrix[b][t]);
    }
  }
  //cout << "ExtNodes: " << nodeVar.size() << endl;
  try{
    //if(show) 
    CPXsetintparam(env,CPX_PARAM_SCRIND,1);
    int status;
    env = CPXopenCPLEX (&status);
    if(status)
      cout << "Error opening environment ... " << status << endl;

    int status_p;
    // Create a cplex network pointer
    prob = CPXNETcreateprob(env, &status_p, "NetFlowProb");
    if(status_p) cout << "Error ... " << status_p << endl;
    // set up a few things for the solver
    // make a call directly to the solver
    // create correct variables first
    int objsen = CPX_MIN;
    CPXData_LR *cd =  getExtNodeExtArcInfo();
    status = CPXNETcopynet(env, prob, objsen, cd->nnodes, &cd->supply[0], NULL, 
			   cd->nedges, &cd->from[0], &cd->to[0], &cd->low[0], &cd->up[0], &cd->obj[0], NULL);
    if(status)
      cout << "Error creating CPXNET ... " << status << endl;
    else
      cout << "\tCreated CPXNET, now optimising ... " << endl;
    delete cd;
  }
  catch(IloException &e){
    cout << endl << "Error during optimisation ..." << e<< endl;
  }
  catch(...){
    cout << endl << "Unknown exception during optimisation ... " << endl;
  }  
}

/*
  Update costs without re-creating the network
  Note, the supplies are scaled: val/1e8
    So the objective needs to be scaled appropriately
*/

int MaxClosure_NetworkFlow_LR::updateSupply(){
  int status;
  CPXData_LR *cd =  getExtNodeExtArcInfo();
  status = CPXNETchgsupply (this->env, this->prob, cd->nnodes, &cd->nodeIds[0], &cd->supply[0]);
  delete cd;
  return status;
}

/*
  Optimise
    Objective must be scaled after
*/

int MaxClosure_NetworkFlow_LR::solve(){
  return CPXNETprimopt(env,prob);
}

/*
  Determine arc and node information from a solution
*/
double MaxClosure_NetworkFlow_LR::getSolInfo(){
  int nnodes = nodeVar.size();
  int nedges = edgeVar.size();
  double obj=0;
  CPXNETgetobjval(env,prob,&obj);
  cout << "\tObjective: " << obj << endl;
  double upper_bound = obj;
  // Get values
  vector<double> x(nedges);
  int status = CPXNETgetx(env, prob, &x[0], 0, nedges-1);
  if(status)
    cout << "Error obtaining arc values ... " << status << endl;

  vector<double> pi(nnodes);
  status = CPXNETgetpi(env, prob, &pi[0], 0, nnodes-1);
  if(status)
    cout << "Error obtaining node values ... " << status << endl;
  // get the values in the solution and set them to the original nodes
  cout << "\tExtNodes:" << endl;
  int count = 0;
  for(int i = 0; i < nnodes; i++){
    if(pi[i] != 0.0) {
      count++;
      setExtNodeVal(i,pi[i]);
    }
  }
  cout << "\tTotal nodes: " << pi.size() << ", nodes used: " << count << endl;
  return upper_bound;
}

/*
  Return the node values after a solve
*/

int MaxClosure_NetworkFlow_LR::getClosure(vector<int>  &sol){
  int nnodes = nodeVar.size();
  vector<double> pi(nnodes);
  int status = CPXNETgetpi(env, prob, &pi[0], 0, nnodes-1);
  if(status)
    cout << "Error obtaining node values ... " << status << endl;
  // get the values in the solution and set them to the original nodes
  //cout << "ExtNodes:" << endl;
  for(int i = 0; i < nnodes; i++){
    if(pi[i] != 0.0) {
      sol.push_back(pi[i]);
    }
  }
  //cout << "Total nodes: " << pi.size() << ", nodes used: " << count << endl;  
  return status;
}

/*
  Set the profit for each node in the network
  Notes:
    The vector "profit" must be the same length as nodeVar
    The vector must also include the profit for the first (dummy) node
*/
void MaxClosure_NetworkFlow_LR::setProfit(const vector<double> &profit){
  for(int i = 0; i < nodeVar.size(); i++){
    nodeVar[i]->setSupply(profit[i]);
  }
}

/*
  Retrieve some information about the nodes
  Note the supplies are scaled. This is required due to precision issues
    with the CPXNET functions
*/

CPXData_LR* MaxClosure_NetworkFlow_LR::getExtNodeExtArcInfo(){
  CPXData_LR *cd = new CPXData_LR();
  int nnodes = nodeVar.size();
  cd->nnodes = nnodes;
  for(int i = 0 ; i < nnodes ; i++ ){
    cd->supply.push_back(nodeVar[i]->getSupply()/1e8);
    //cd->supply.push_back(nodeVar[i]->getSupply());
    cd->nodeIds.push_back(nodeVar[i]->getID());
  }
  int nedges = edgeVar.size();
  cd->nedges = nedges;
  for(int i = 0 ; i < nedges ; i++ ){
    cd->from.push_back(edgeVar[i]->getSrc()->getID());
    cd->to.push_back(edgeVar[i]->getTgt()->getID());
    cd->low.push_back(edgeVar[i]->getLb());
    cd->up.push_back(edgeVar[i]->getUb());
    cd->obj.push_back(edgeVar[i]->getCost());
    cd->edgeIds.push_back(edgeVar[i]->getID());
  }
  return cd;
}

/*
  Set values for the nodes
 */

void MaxClosure_NetworkFlow_LR::setExtNodeVal(int i, double value){
  int blk = nodeIdMapping[i].first;
  int time = nodeIdMapping[i].second;
  NetExtNode *nd = nodeVarMatrix[blk][time];
  nd->setValue(value);
}

/*
  Get a node's (MILP determined) value
*/
double MaxClosure_NetworkFlow_LR::getExtNodeVal(int block, int time){
  return nodeVarMatrix[block][time]->getValue();
}

/*
  Get the supply values
*/


void MaxClosure_NetworkFlow_LR::getSupplyInfo(pair<vector<int>, vector<double> > &supply){
  vector<double> supp;
  vector<int> ids;
  //supp.push_back(dummy->getSupply());
  int nnodes = nodeVar.size();
  for(int i = 0 ; i < nnodes ; i++ ){
    supp.push_back(nodeVar[i]->getSupply());
    ids.push_back(nodeVar[i]->getID());
  }
  supply.first = ids;
  supply.second = supp;
}

/*
  Reset values after a solve
*/

void MaxClosure_NetworkFlow_LR::resetValues(){
  for(int i = 0; i < nodeVar.size() ; i++){
    nodeVar[i]->setValue(0.0);
  }
}


void MaxClosure_NetworkFlow_LR::updateProfit(const vector<vector<double> > &profit){

  this->profit.clear();
  this->profit = profit;

  setProfitSum();
  // The actual update
  dummy->setSupply(profitSum);
  for(int b = 0; b < nodeVarMatrix.size(); b++ ){
    for(int t = 0; t < nodeVarMatrix[b].size(); t++ ){
      nodeVarMatrix[b][t]->setSupply(-profit[b][t]);
    }
  }  
}

void MaxClosure_NetworkFlow_LR::setProfitSum(){
  this->profitSum = 0.0;
  int nB = profit.size();
  int t_max = profit[0].size();
  for(int b = 0 ; b < nB ; b++) {
    for(int t=0; t<t_max; t++){
      this->profitSum += this->profit[b][t];
    }
  }
}

/*
  Clean up the memory allocated
*/
MaxClosure_NetworkFlow_LR::~MaxClosure_NetworkFlow_LR(){

}




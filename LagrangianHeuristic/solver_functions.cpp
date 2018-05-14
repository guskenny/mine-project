/***************************************************************************
                            SolverFunctions.cpp
                         ----------------------------------
    last modified   : 13/5/2008
    copyright       : (C) 2013 by Dhananjay Thiruvady
    libraries	    : .
    description	    : Implementation of functions for Solver
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 
#include <iostream>
#include <fstream>
#include <list>
#include <cmath>

using namespace std;

#include "solver_functions.h"

/*
  Construct the network flow class
*/

SolverFunctions::SolverFunctions(Daten &data, int diff_amt){
  this->data = &data;
  this->diff_amt = diff_amt;

  setMaxBlockProfit();
  setSucc();
  setBackPropProfit();
  setNegBlock();
  //setDeepBackPropProfit();
  //setBackPropProfitAvg();
  //setPosBackPropProfit();
}

/*
  Update the Lagrangian multipliers
  Currently, this is a straightforward subgradient optimisation
  This is to be improved by the volume algorithm or bundle methods
*/

vector<vector<double> > SolverFunctions::update_lambda(double gap, const vector<vector<double> > &lambda, double delta, const vector<vector<vector<double> > > &Yvar){
  // set up some parameters
  const char pT = data->getProbType();
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const int r_max = data->getnResources();
  const double rate = data->getDiscountRate();
  std::vector<Block> * blocks=data->getBlock();

  vector<vector<double> > updated_lambda;
  double delta_kt = 0.0;
  for(int r=0;r<r_max;r++){
    for(int t=0;t<t_max;t++){
      double factor=0.0;
      for(int b=0; b<nB; b++){
	for(int d=0; d<d_max; d++){
	  double coef = (*blocks)[b].getRCoef(d,r);
	  factor += coef*Yvar[b][t][d];
	}
      }
      // the subgradient depends on the type of constraint we are working with
      char cType = data->getResConstrType(r, t);
      if(cType == 'L'){ 
	factor-=data->getLimit(r, t);
      }
      else if(cType == 'R'){ // we are not dealing with this case for now
	// factor+=data->getLimit(r, t);
      }
      else if(cType == 'I'){// we are not dealing with this case for now
	//factor+=data->getLimit(r, t, 1);
	//factor-=data->getLimit(r, t);
      }
      delta_kt += factor*factor;
    }
  }
  updated_lambda = lambda;
  for(int r=0;r<r_max;r++){
    for(int t=0;t<t_max;t++){
      double factor=0.0;
      for(int b=0; b<nB; b++){
	for(int d=0; d<d_max; d++){
	  double coef = (*blocks)[b].getRCoef(d,r);
	  factor += coef*Yvar[b][t][d];
	}
      }
      // the subgradient depends on the type of constraint we are working with
      char cType = data->getResConstrType(r, t);
      if(cType == 'L'){ 
	factor-=data->getLimit(r, t);
      }
      else if(cType == 'R'){ // we are not dealing with this case for now
	//factor+=data->getLimit(r, t);
      }
      else if(cType == 'I'){ // we are not dealing with this case for now
	//factor+=data->getLimit(r, t, 1);
	//factor-=data->getLimit(r, t);
      }
      updated_lambda[r][t] += (delta*gap*factor)/delta_kt;
      if(updated_lambda[r][t]<0.0) updated_lambda[r][t] = 0.0;
      // validation purposes
      //cout << r << "\t" << t << "\t" << updated_lambda[r][t] << endl;
    }
  }
  return updated_lambda;
} 

/*
  Compute an average block profit
  // The following description is from Andreas:
  //  "Given a multiplier, calculate lagrangian vector as follows:
  //   avgBlockProfit = avg( max(profitBase[b][d] for d) for b)
  //   lag[r][0] = mult * avgBlockProfit / avg(q[b][d][r] for d for b) / len(R)
  //   lag[r][t] = lag[r][0] / (1+self.discount)**t
  //   This should mean that at mult=1 average reduced profit is about 0
*/

double SolverFunctions::ComputeAvgBlockProfitSum(){  
  // initialise some parameters from the data
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  std::vector<Block> * blocks=data->getBlock();

  // Find the average block profit
  double BlockProfitSum = 0.0;
  for(int b=0; b<nB; b++){
    // first determine the desination with the max profit
    double max_profit = -1e99;
    for(int d=0; d<d_max; d++){
      double t_profit = (*blocks)[b].getProfit(d); 
      if(t_profit > max_profit) {
	max_profit = t_profit;
      }
    }
    BlockProfitSum += max_profit;
  }
  return BlockProfitSum/nB;
}

/*
  Compute average resource use
*/

double SolverFunctions::ComputeAvgResourceUse(int r){  
  // initialise some parameters from the data
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  std::vector<Block> * blocks=data->getBlock();

  // Find the average block profit
  double ResourceSum = 0.0;
  for(int b=0; b<nB; b++){
    // first determine the desination with the max profit
    double max_resource = -1e99;
    for(int d=0; d<d_max; d++){
      double t_resource = (*blocks)[b].getRCoef(d,r);
      if(t_resource > max_resource) {
	max_resource = t_resource;
      }
    }
    ResourceSum += max_resource;
  }
  //cout << "Resource sum:" << ResourceSum << endl;
  return ResourceSum/nB;
}


/*
  Network related functions
 */

void SolverFunctions::setProfit(   vector<double> &profit,
			 vector<vector<double> > &profit_time){

  const int nB = data->getNBlock();
  const int d_max = data->getnDestination();
  const int t_max = data->getNPeriod();
  std::vector<Block> * blocks=data->getBlock();

  profit.resize(nB,0.0);
  for(int b=0; b<nB; b++){
    // determine the desination with the max profit
    double max_profit = -1e99;
    for(int d=0; d<d_max; d++){
      double t_profit = (*blocks)[b].getProfit(d); 
      if(t_profit > max_profit) {
	max_profit = t_profit;
      }
    }
    profit[b] = max_profit;
  }
  // set up a dummy profit time, needed for other structures
  profit_time.resize(nB);
  for(int b=0; b<nB; b++){
    profit_time[b].resize(t_max,0.0);
  }
}

/*
  Set up the profit vector considering the Lagrangian multipliers
  Note, use this function only if using the Lagrangian heuristic
*/

void SolverFunctions::setLambdaProfit(  vector<vector<double> > &profit_time,
			       const vector<vector<double> > &lambda){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const int r_max = data->getnResources();
  const double rate = data->getDiscountRate();
  std::vector<Block> * blocks=data->getBlock();

  // Set up profit with a time component
  profit_time.resize(nB);
  for(int b=0; b<nB; b++){
    profit_time[b].resize(t_max,0.0);
    for(int t=0; t<t_max; t++){
      // first determine the desination with the max profit
      double max_profit = -1e99;
      for(int d=0; d<d_max; d++){
	double t_profit = (*blocks)[b].getProfit(d)/pow(1+rate,t); 
	for(int r=0; r<r_max; r++)
	  t_profit -= (lambda[r][t] * (*blocks)[b].getRCoef(d,r));
	if(t_profit > max_profit) {
	   max_profit = t_profit;
	}
      }
      profit_time[b][t] = max_profit;
    }
  }
  // adjust profits
  for(int b = 0 ; b < nB ; b++) {
    for(int t=0; t<t_max-1; t++){
      profit_time[b][t] -= profit_time[b][t+1];
    }
  }
}

/*
  Set max block profit
*/
void SolverFunctions::setMaxBlockProfit(){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const int r_max = data->getnResources();
  const double rate = data->getDiscountRate();
  std::vector<Block> * blocks=data->getBlock();

  this->max_block_profit.resize(nB);
  for(int b=0; b<nB; b++){
    this->max_block_profit[b].resize(t_max,0);
    for(int t=0; t<t_max; t++){
      // first determine the desination with the max profit
      double max_profit = -1e99;
      int dest = 0;
      for(int d=0; d<d_max; d++){
	double t_profit = (*blocks)[b].getProfit(d)/pow(1+rate,t); 
	if(t_profit > max_profit) {
	   max_profit = t_profit;
	   dest = d;
	}
      }
      this->max_block_profit[b][t] = dest;
    }
  }
}

/*
  Generate a solution like Yvar
*/

vector<vector<vector<double> > > SolverFunctions::generateSolution(MaxClosure_NetworkFlow_LR *net){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();

  vector<vector<vector<double> > > sol(nB);
  for(int b=0; b<nB; b++){
    sol[b].resize(t_max);
    for(int t=0; t<t_max; t++){
      sol[b][t].resize(d_max,0.0);
    }
  }

  // Generate a solution
  for(int b = 0; b < nB; b++){
    for(int t = 0; t < t_max; t++){
      if(net->getExtNodeVal(b,t) != 0.0){
	int dest = max_block_profit[b][t];
	sol[b][t][dest] = 1.0;
	break;
      }
    }
  }
  return sol;
}

/*
  Initialise lambda
*/

void SolverFunctions::setLambda(vector<vector<double> > &lambda){
  const int t_max = data->getNPeriod();
  const int r_max = data->getnResources();
  const double rate = data->getDiscountRate();

  // Initial multipliers, set them up like Andreas does
  double mult = 0.1;
  double avgBlockProfit = ComputeAvgBlockProfitSum();
  //setMaxBlockProfit(data);
  cout << "\tAverage profit: " << avgBlockProfit << endl;
  double min_lag = 1e99;
  double max_lag = -1e99;

  lambda.resize(r_max);
  for(int r = 0; r < r_max ; r++){
    lambda[r].resize(t_max,0.0);
    lambda[r][0] = mult * avgBlockProfit / ComputeAvgResourceUse(r) / r_max;
    if(lambda[r][0] > max_lag) max_lag = lambda[r][0];
    if(lambda[r][0] < min_lag) min_lag = lambda[r][0];
    for(int t = 1; t < t_max ; t++){
      lambda[r][t] = lambda[r][0] / pow(2 + rate,t);   
      if(lambda[r][t] > max_lag) max_lag = lambda[r][t];
      if(lambda[r][t] < min_lag) min_lag = lambda[r][t];
    }
  }
  cout << "\tMultiplier range: (" << min_lag << "," << max_lag << ")" << endl; 
}

/*
  Compare solutions
 */

bool SolverFunctions::validateSols(const vector<vector<vector<double> > > &Yvar,
				   const vector<vector<vector<double> > > &Yvar_mod){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  bool same = true;
  // check out what has been set in Yvar
  for(int i =0;i<nB;i++){
    for(int j =0;j<t_max;j++){
      for(int k =0;k<d_max;k++){
	if(Yvar[i][j][k] > 0.01){
	  //cout << "Block " << i << ", time " << j << ", dest " << k << endl;
	}
	if(Yvar[i][j][k] != Yvar_mod[i][j][k]){
	  same = false;
	  cout << "Initial problem ... " << i << ":" << j << ":" << k<< endl;
	}
      }
    }
  }
  return same;
}


/*
  Set the profit per block
*/

void SolverFunctions::setProfit(){
  const int nB = data->getNBlock();
  const int d_max = data->getnDestination();
  const int t_max = data->getNPeriod();
  std::vector<Block> * blocks=data->getBlock();
  const int r_max = data->getnResources();

  profit.resize(nB,0.0);
  block_dest.resize(nB,0);
  for(int b=0; b<nB; b++){
    // determine the desination with the max profit
    double max_profit = -1e99;
    int dest = 0;

    // first compute resources for each destination
    vector<double> ru(d_max,0.0);
    for(int d=0; d<d_max; d++){
      for(int r = 0; r < r_max; r++){
	double coef = (*blocks)[b].getRCoef(d,r);
	ru[d] += coef;
      }
    }
    for(int d=0; d<d_max; d++){
      double t_profit = (*blocks)[b].getProfit(d); 
      if(abs(t_profit-max_profit) < diff_amt){ // nearly equal profits, choose the less resource intensive destination
	double min_resource = 1e99;
	int dest1 = 0;
	for(int d1 = 0; d1 < d_max; d1++){
	  if(ru[d1] < min_resource){
	    min_resource = ru[d1];
	    dest1 = d1;
	  }
	}	
	max_profit = (*blocks)[b].getProfit(dest1); ;
	dest = dest1;
      }
      else if(t_profit > max_profit) {
	max_profit = t_profit;
	dest = d;
      }
    }
    profit[b] = max_profit;
    block_dest[b] = dest;
  }
}

/*
  Set a profit based on back propagation
*/

void SolverFunctions::setBackPropProfit(){
  const int nB = data->getNBlock();
  const int r_max = data->getnResources();
  std::vector<Block> * blocks=data->getBlock();
  const double rate = data->getDiscountRate();
  setProfit();
  vector<double> t_prof = this->profit;
  //vector<double> tt_prof = this->profit;
  // now update the profits 
  // go over the successors
  for(int blk = 0; blk < nB; blk++){
    for (set<int>::iterator lit = succ[blk].begin(); lit != succ[blk].end(); lit++) { 
      int t_blk = (*lit);
      int n = (*blocks)[t_blk].getNumPred();
      if(t_prof[t_blk] < 0.0 && succ[t_blk].size() == 0 )  // this block will be discarded anyway
	  continue;
      //if(t_prof[t_blk] < 0.0)  // this block will be discarded anyway
	  //  continue;
      // get resource usage for that destination
      profit[blk] += t_prof[t_blk]/(n+1);
    }
    /*double r_cons = 1.0;
    for(int r=0; r<r_max; r++){
      double coef = (*blocks)[blk].getRCoef(block_dest[blk],r);
      r_cons += coef;
    }
    profit[blk] += r_cons;*/
  }
  /*for(int blk = 0; blk < nB; blk++){
    for (set<int>::iterator lit = succ[blk].begin(); lit != succ[blk].end(); lit++) { 
      int t_blk = (*lit);
      profit[blk] += tt_prof[t_blk];
    }
    }*/

}

/*
  Set a profit based on back propagation
*/

void SolverFunctions::setPosBackPropProfit(){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const double rate = data->getDiscountRate();
  setProfit();
  vector<double> t_prof = this->profit;
  // now update the profits 
  // go over the successors
  for(int blk = 0; blk < nB; blk++){
    for (set<int>::iterator lit = succ[blk].begin(); lit != succ[blk].end(); lit++) { 
      int t_blk = (*lit);
      if(t_prof[t_blk] < 0.0)
	continue;
      //profit[blk] += t_prof[t_blk]/pow(1+rate,t_max-1);
      //else
      profit[blk] += t_prof[t_blk];
    }
  }
}

/*
  Set a profit based on back propagation
  The average is across time
*/

void SolverFunctions::setBackPropProfitAvg(){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const double rate = data->getDiscountRate();
  setProfit();
  vector<double> t_prof = this->profit;
  // now update the profits 
  // go over the successors
  for(int t = 0; t < t_max; t++){
    for(int blk = 0; blk < nB; blk++){
      for (set<int>::iterator lit = succ[blk].begin(); lit != succ[blk].end(); lit++) { 
	int t_blk = (*lit);
	profit[blk] += t_prof[t_blk]/pow(1+rate,t);
      }
    }
  }
  // normalise
  for(int blk = 0; blk < nB; blk++){
    profit[blk] /= t_max;
  }
}

/*
  Recursive function to get successor profits
  No need for a base case as it returns 0 if there are no successors
  The second arguement (vector) is used for memoization
 */
double SolverFunctions::getSuccProfit(int block, vector<pair<double, bool> > &mem_profit){
  const int nB = data->getNBlock();
  std::vector<Block> * blocks=data->getBlock();
  const double rate = data->getDiscountRate();
  double prof = 0.0;
  for (set<int>::iterator lit = succ[block].begin(); lit != succ[block].end(); lit++) { 
    //mem_profit[(*lit)].first += profit[block]/nB;;
    int n = (*blocks)[(*lit)].getNumPred();
    if(mem_profit[(*lit)].second == false)
      prof += getSuccProfit((*lit), mem_profit);///(n+1)/rate;
    //prof += max(0.0,getSuccProfit((*lit), mem_profit)/(n+1)/rate);
    else 
      prof += mem_profit[(*lit)].first;///(n+1)/rate;
    //prof += max(0.0,mem_profit[(*lit)].first/(n+1)/rate);
  }  
  double factor = prof;
  if(prof < 0.0) factor /= (nB*rate);
  if(profit[block] < 0.0 && succ[block].size()==0)
      mem_profit[block].first += 0.0;  
  else
    mem_profit[block].first +=  factor + profit[block];  
  mem_profit[block].second = true;  
  return mem_profit[block].first;
}

/*
  Set a profit based on back propagation to the bottom of the pit
  Propagate only profit/n where n are the immediate predecessors
*/

void SolverFunctions::setDeepBackPropProfit(){
  const int nB = data->getNBlock();
  const int r_max = data->getnResources();
  std::vector<Block> * blocks=data->getBlock();
  setProfit();
  vector<double> t_prof(nB,0.0);
  vector<pair<double,bool> > mem_profit(nB,make_pair(0.0,false));
  // now update the profits 
  // go over the successors
  for(int blk = 0; blk < nB; blk++){
    if(mem_profit[blk].second == false)
      t_prof[blk] += getSuccProfit(blk, mem_profit);
    //t_prof[blk] += max(0.0,getSuccProfit(blk, mem_profit));
    else
      t_prof[blk] += mem_profit[blk].first;
    //t_prof[blk] += max(0.0,mem_profit[blk].first);
    //cout << blk << ", " << t_prof[blk] << ", " << profit[blk] << endl; 
    /*double r_cons = 1.0;
    for(int r=0; r<r_max; r++){
      double coef = (*blocks)[blk].getRCoef(block_dest[blk],r);
      r_cons += coef;
    }
    t_prof[blk] += r_cons;*/
  }
  //for(int blk = 0; blk < nB; blk++){
    //if(t_prof[blk]-profit[blk] > 0)
      //cout << blk << ":" << t_prof[blk]-profit[blk] << endl;
  //}
  this->profit.clear();
  this->profit = t_prof;
}

/*void SolverFunctions::setDeepBackPropProfit(){
  const int nB = data->getNBlock();
  std::vector<Block> * blocks=data->getBlock();
  const double rate = data->getDiscountRate();
  setProfit();
  vector<double> t_prof = this->profit;
  vector<int> succ_size(nB,0);
  list<int> avail_blocks;
  for(int blk = 0; blk < nB; blk++){
    succ_size[blk] += succ[blk].size();
    std::vector<int> * pred= (*blocks)[blk].getPreds();
    if((*pred).size() == 0){
      avail_blocks.push_back(blk);
    }
  }
  // until all blocks are finished
  //int level = 0;
  while(avail_blocks.size() > 0){
    list<int> poss_blocks;
    for (list<int>::iterator lit = avail_blocks.begin(); lit != avail_blocks.end(); lit++) { 
      //cout << (*lit) << ":" << level << endl;
      std::vector<int> * pred= (*blocks)[(*lit)].getPreds();
      for(int p=0; p<(*pred).size(); p++){
	int b = (*pred)[p];
	if(t_prof[(*lit)] > 0.0)
	  t_prof[b] += t_prof[(*lit)]/(*pred).size()/rate;
	succ_size[b]--;
	if(succ_size[b] == 0){
	  //cout << b << " is available " << endl;
	  poss_blocks.push_back(b);
	}
      }
    }
    avail_blocks.clear();
    avail_blocks = poss_blocks;
    //level++;
  }
  for(int blk = 0; blk < nB; blk++){
    if(t_prof[blk]-profit[blk] > 0)
      cout << blk << ":" << t_prof[blk]-profit[blk] << endl;
  }
  this->profit.clear();
  this->profit = t_prof;
  }*/

/*
  Determine the successors for each block
*/

void SolverFunctions::setSucc(){
  const int nB = data->getNBlock();
  std::vector<Block> * blocks=data->getBlock();
  // determine the successors
  succ.resize(nB);
  for(int blk = 0; blk < nB; blk++){
    std::vector<int> * pred= (*blocks)[blk].getPreds();
    for(int p=0; p<(*pred).size(); p++){
      int b = (*pred)[p];
      succ[b].insert(blk);
    }  
  }
}

/*
  Recursive function to determine negative successors 
*/

bool SolverFunctions::getNegSucc(int block, vector<pair<bool, bool> > &mem_succ){
  bool ns = true;
  if(profit[block] > 0.0) ns = false;
  else{
    for (set<int>::iterator lit = succ[block].begin(); lit != succ[block].end() && ns; lit++) { 
      if(mem_succ[(*lit)].second == false)
	ns = getNegSucc((*lit), mem_succ);
      else 
	ns = mem_succ[(*lit)].first;
    }  
  }
  mem_succ[block].first = ns;
  mem_succ[block].second = true;  
  return mem_succ[block].first;
}

/*
  Recursive function to determine negative chains
  Note: this is only valid at time point 0
*/

double SolverFunctions::getNegSucc(int block, vector<pair<double, bool> > &mem_succ){
  if(mem_succ[block].second == false){ // we don't have a value for the block
    double prof = 0.0;
    for (set<int>::iterator lit = succ[block].begin(); lit != succ[block].end(); lit++) { 
      if(mem_succ[(*lit)].second == false)
	prof += getNegSucc((*lit), mem_succ);
      else 
	prof += mem_succ[(*lit)].first;
    }  
    mem_succ[block].first = prof;
    mem_succ[block].second = true;  
  }
  return mem_succ[block].first;
}

/*
  Display some information
*/
void SolverFunctions::dispSucc(int block, vector<bool> &visited){
  if(visited[block]) 
    return;
  cout << "\tBlock: " << block << ", profit: " << profit[block] << endl;
  for (set<int>::iterator lit = succ[block].begin(); lit != succ[block].end(); lit++) { 
    dispSucc((*lit),visited);
  }
  visited[block] = true;
}

/*
  Set the negative blocks
*/

void SolverFunctions::setNegBlock(){
  const int nB = data->getNBlock();
  std::vector<Block> * blocks=data->getBlock();
  // determine the negative blocks
  negBlocks.resize(nB,false);
  vector<pair<bool,bool> > mem_succ(nB,make_pair(false,false));
  // go over the successors
  for(int blk = 0; blk < nB; blk++){
    bool t_succ = false;
    if(mem_succ[blk].second == false)
      t_succ = getNegSucc(blk, mem_succ);
    else
      t_succ = mem_succ[blk].first;
    negBlocks[blk] = t_succ;
  }

  // Remove chains
  /*vector<pair<double,bool> > mem_succ(nB,make_pair(-1e30,false));
  // go over the successors
  for(int blk = 0; blk < nB; blk++){
    double t_succ = false;
    if(mem_succ[blk].second == false)
      t_succ = getNegSucc(blk, mem_succ);
    else
      t_succ = mem_succ[blk].first;
    if(mem_succ[blk].first < 0.0)
      negBlocks[blk] = true;
      }*/

  // Validate the negative paths
  /*for(int blk = 0; blk < nB; blk++){
    if(negBlocks[blk] == true){
      cout << "Block: " << blk << " has a negative path with profit: " << profit[blk] << endl;
      vector<bool> visited(nB, false);
      dispSucc(blk,visited);
    }
  }
  cout << "Hard test: " << endl;
  vector<bool> visited(nB, false);
  dispSucc(8540, visited);*/
}

/*
  Generate a solution needed by the parallel side of things
  Note the objective is not needed here
 */
void SolverFunctions::generateParSolution(ACO_Solution *a, Sol_Int &sol){
  const int t_max = data->getNPeriod();
  // Do the conversions
  const int d_max = data->getnDestination();  
  vector<int> times = a->getBlockTimes();
  vector<int> dests = a->getBlockDests(); // this is limiting in the case of the MIP
  vector<int> X;
  vector<vector<double> > Y;
  
  int nB = times.size();
  X.clear();
  Y.clear();
  X.resize(nB,t_max);
  Y.resize(nB);
  
  for(int b=0;b<nB;b++){
    if(times[b] != -1) 
      X[b] = times[b];
    Y[b].resize(d_max,0.0);
    int dest = dests[b];
    Y[b][dest] = 1.0;
  }
  sol.x = X;
  sol.y = Y;
  sol.obj = a->getObj();
  sol.nT = t_max;
}

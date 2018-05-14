/*************************************************************************** 
                            solution.cpp 
                         ------------------- 
    last modified   : 13/5/2016 
    copyright       : (C) 2016 by Dhananjay Thiruvady 
    libraries		: . 
    description		: a solution to the probelm 
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
#include <new> 
#include <vector> 
#include <cmath> 
#include <list>
 
using namespace std; 
 
#include "ACO_solution.h" 
#include "pheromones.h" 
#include "Timer.h"
 
ACO_Solution::ACO_Solution(Daten *data, const vector<double> &profit, 
			   const vector<int> &bd, const vector<set<int> > &succ){ 
  const int nB = data->getNBlock();
  /*const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  */
  cost=0.0; 
  violations=0; 
  block_times.resize(nB,-1);
  block_seq.resize(nB,-1);
  this->block_dest = bd;
  this->data = data;
  this->profit = profit;
  this->succ = succ;
  /*this->Y.resize(nB);
  for(int b = 0; b < nB; b++){
    this->Y[b].resize(t_max);
    for(int t = 0; t < t_max; t++){
      this->Y[b][t].resize(d_max,0.0);
    }
    }*/
}

/*
  A second constructor to map a problem data solution into an ACO solution
  Also sets an initial starting state using he heuristic
  Note: bd or default block destinations are not used since Y should have 
     the ideal destinations for each block. However, it is passed anyway as 
     there maybe some reason for using it at a later stage
*/

ACO_Solution::ACO_Solution(Daten *data, const vector<double> &profit, 
			   const vector<int> &bd, const vector<set<int> > &succ,
			   const vector<bool> &negBlocks,
			   const vector<int> &X,
			   const vector<vector<double> > &Y,
			   bool exact
			   ){ 
  const int nB = data->getNBlock();
  cost=0.0; 
  violations=0; 
  block_times = X;
  //block_times.resize(nB, 0);
  block_seq.resize(nB,-1);
  // Should be set from Y in final version
  this->block_dest = bd;
  this->data = data;
  this->profit = profit;
  this->succ = succ;
  if(exact) this->setBlocksExact();
  else this->setBlocksWithPrec(negBlocks);
}

ACO_Solution* ACO_Solution::copy(){
  ACO_Solution *sol = new ACO_Solution(this->data,this->profit,this->block_dest,this->succ);
  sol->cost = this->cost;
  sol->violations = this->violations;
  sol->block_times = this->block_times;
  sol->block_seq = this->block_seq;
  return sol;
}

/*
  Network related functions
 */


/*
  Check if the preceding jobs have been completed
*/

bool ACO_Solution::checkPrecedence(int blk, const vector<bool> &blocks_completed){ 
  std::vector<Block> * blocks=data->getBlock();
  std::vector<int> * pred = (*blocks)[blk].getPreds();
  bool not_done = false;
  for(int p=0; p<(*pred).size() && !not_done; p++){
    int b = (*pred)[p];
    if(!blocks_completed[b]) not_done = true;
  }
  return not_done;
} 

/*
  Determine times from the Yvar solution
*/
void ACO_Solution::setTimes(const vector<vector<vector<double> > > &Yvar){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();

  // Set destinations
  // Currently (6/7/16), blocks can only be sent to one destination
  // This should change to proportions, but needs some thinking
  for(int b=0 ; b < nB ; b++){
    bool dest_done = false;
    for(int d=0 ; d < d_max && !dest_done; d++){
      for(int t=0 ; t < t_max && !dest_done; t++){
	if(Yvar[b][t][d] > 0){
	  block_dest[b] = d;
	  dest_done = true;
	}
      }
    }
  }

  for(int b=0 ; b < nB ; b++){
    bool time_done = false;
    for(int t=0 ; t < t_max && !time_done; t++){
      for(int d=0 ; d < d_max && !time_done; d++){
	if(Yvar[b][t][d] > 0){
	  block_times[b] = t;
	  time_done = true;
	}
      }
    }
  }
  setSequence(Yvar);
}

/*
  Determine times from the Yvar solution
*/
void ACO_Solution::setTimes(const vector<int> &times){
  this->block_times.clear();
  this->block_times = times;
}


/*
  Sort a sequence
  The profit vector needs to be modified
*/

void ACO_Solution::sortSeq(vector<vector<int> > &timeBlocks){
  // set up some parameters
  const int t_max = data->getNPeriod();
  vector<double> profits = this->profit;
  // At each time point we now have all the correct blocks
  // Set their ordering according profit values first
  for(int t = 0; t < t_max; t++ ){
    if(timeBlocks[t].size() < 2) // nothing to sort
      continue;
    // sort the blocks in order of profit
    for(int b = 0; b < timeBlocks[t].size()-1; b++){
      int blk = timeBlocks[t][b];
      double maxProfit = profits[blk];
      int index = b;
      for(int b1 = b+1; b1 < timeBlocks[t].size(); b1++){
	int blk1 = timeBlocks[t][b1];
	if(maxProfit < profits[blk1]){
	  maxProfit = profits[blk1];
	  index = b1;
	}
      }
      // Swap index with b
      int temp = timeBlocks[t][b];
      timeBlocks[t][b] = timeBlocks[t][index];
      timeBlocks[t][index] = temp; 
    }
  }

  // Set the sequence
  int bs_index = 0;
  for(int t = 0; t < t_max; t++ ){
    //cout << "Time: " << t << endl;
    for(int b = 0; b < timeBlocks[t].size(); b++){
      int blk = timeBlocks[t][b];
      this->block_seq[bs_index] = blk;
      bs_index++;
      //cout << "Block: " << blk << ", profit: " << profits[blk] << " index: " << bs_index << endl;
    }
  }
}

/*
  Sort a sequence probabilistically
  I.e., assign the highest probability to the highest value job, etc
  Do a roulette wheel selection based on the above
  The profit vector needs to be modified, has a bug but doesn't work well
*/

void ACO_Solution::sortProbSeq(vector<vector<int> > &timeBlocks){
  // set up some parameters
  const int t_max = data->getNPeriod();
  vector<double> profits = this->profit;
  // determine the minimum profit
  double min_prof = 1e30;
  for(int b = 0; b < profits.size(); b++){
    if(profits[b] < min_prof) min_prof = profits[b];
  }
  
  time_t t;
  Random *rnd;
  rnd = new Random((unsigned) time(&t));
  rnd->next();
  // At each time point we now have all the correct blocks
  // Set their ordering according profit values first
  for(int t = 0; t < t_max; t++ ){
    if(timeBlocks[t].size() < 2) // nothing to sort
      continue;
    // sort the blocks in order of profit
    vector<int> seq = timeBlocks[t];
    vector<bool> done(seq.size(),false);
    for(int b = 0; b < seq.size(); b++){
      vector<double> selProb(seq.size(),0.0);
      double sumProbs = 0.0;
      for(int b1 = 0; b1 < seq.size(); b1++){
	int blk = seq[b1];	
        if(done[b1]) continue;
	selProb[b1] += profits[blk]+min_prof;
	sumProbs += selProb[b1]; 
      }
      double r_num= rnd->next();   
      double r = sumProbs*r_num;	 
      int block = 0; 
      double prob = selProb[block]; 
      while( prob < r ){ 
	block++; 
	prob += selProb[block]; 
      }	 
      done[block] = true;
      timeBlocks[t][b] = seq[block];
    }
  } 
  delete rnd;
  // Set the sequence
  int bs_index = 0;
  for(int t = 0; t < t_max; t++ ){
    //cout << "Time: " << t << endl;
    for(int b = 0; b < timeBlocks[t].size(); b++){
      int blk = timeBlocks[t][b];
      this->block_seq[bs_index] = blk;
      bs_index++;
      //cout << "Block: " << blk << ", profit: " << profits[blk] << " index: " << bs_index << endl;
    }
  }
}

/* 
   Determine a sequence sorted by time and max profit
*/

void ACO_Solution::setSequence(const vector<vector<vector<double> > > &Yvar){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();

  // Object which holds which blocks go at which time
  vector<vector<int> > timeBlocks(t_max);
  for(int t = 0; t < t_max; t++ ){
    for(int b = 0; b < nB; b++){
      for(int d=0 ; d < d_max; d++){
      	if(Yvar[b][t][d] > 0.0){
	  timeBlocks[t].push_back(b);
	  break;
	}
      }
    }    
  }
  sortSeq(timeBlocks);
}

/* 
   Determine a sequence sorted by time and max profit
   X is a 0-1 vector of block by time where X_bt has a value of 1
     if the block b is mined at time t 
*/

void ACO_Solution::setSequence(const vector<vector<double> > &X){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();

  // Object which holds which blocks go at which time
  vector<vector<int> > timeBlocks(t_max);
  for(int t = 0; t < t_max; t++ ){
    for(int b = 0; b < nB; b++){
      if(X[b][t] > 0.0){
	timeBlocks[t].push_back(b);
      }
    }    
  }
  sortSeq(timeBlocks);
  //sortProbSeq(timeBlocks);
}

/*
  Assumes preceding jobs have been completed
  Returns the earliest time a block can be scheduled given predecessors
  Shift is used for negative valued blocks
*/

int ACO_Solution::earliestTime(int b, int shift){ 
  std::vector<Block> * blocks=data->getBlock();
  std::vector<int> * pred = (*blocks)[b].getPreds();
  int n = (*blocks)[b].getNumPred();
  int earliest = 0;
  for(int p=0; p<n; p++){
    int b = (*pred)[p];
    if(block_times[b] > earliest) earliest = block_times[b];
  }
  //if(profit[b] < 0.0 && earliest < t_max/2) earliest = t_max/2;
  if(profit[b] < 0.0) earliest += shift;
  return earliest;
} 

/*
  Set the blocks here
  Also determine the objective
*/

double ACO_Solution::setBlocks(const vector<bool> &negBlocks, int shift){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const int r_max = data->getnResources();
  const double rate = data->getDiscountRate();
  std::vector<Block> * blocks=data->getBlock();

  double objective = 0.0;
  // now set the blocks with precedences
  vector<int> blocks_waiting;
  //vector<bool> blocks_completed(nB,false);

  // set an artificial horizon
  // note, the objective will only include blocks within the original horizon
  int horizon = t_max;
  vector<vector<double> > resources_used(horizon);

  for(int t = 0; t < horizon; t++) {
    resources_used[t].resize(r_max,0.0);
  }
  // determine precedence vector
  vector<int> prec_vector(nB,0);
  for(int b = 0; b < nB; b++){
    prec_vector[b] = (*blocks)[b].getNumPred();    
  }

  vector<vector<vector<double> > > Y(nB);
  for(int b = 0; b < nB; b++){
    Y[b].resize(t_max);
    for(int t = 0; t < t_max; t++){
      Y[b][t].resize(d_max,0.0);
    }
  }
  int c1 = 0;
  vector<bool> bd(nB,false);

  for(int b = 0; b < block_times.size(); b++){
    int blk = block_seq[b]; // this is needed as blk will change to something on the waiting list
    if(blk == -1) {
      continue;
    }
    // check for successors and negative value
    // also check for negaitve paths
    if((succ[blk].size() == 0 && profit[blk] < 0) || negBlocks[blk]) {
    //if((succ[blk].size() == 0 && profit[blk] < 0)) {
      //block_seq[b] = -1;
      block_times[blk] = -1;
      //cout << "Useless block: " << blk << " with profit: " << profit[blk] << endl;
      continue;
    }
    //bd[blk] = true;
    //bool cant_schedule=checkPrecedence(blk,blocks_completed); 
    if(prec_vector[blk] > 0) {
    //if(cant_schedule){
      blocks_waiting.push_back(blk);
    }
    else{ 
      bool blocks_found = true; 
      while(blocks_found){ 
	// update prec_vector for all successors of blk
	for (set<int>::iterator lit = succ[blk].begin(); lit != succ[blk].end(); lit++) { 
	  int t_blk = (*lit);
	  prec_vector[t_blk]--;
	}

	//blocks_completed[blk]=true; 
	bool valid_block = false; 
	int t = earliestTime(blk, shift); 
	while(t < horizon && !valid_block){
	  bool block_ok = true; 
	  // assuming CType = 'L'
	  for(int r = 0; r < r_max; r++){
	    int d = block_dest[blk];
	    double total = (*blocks)[blk].getRCoef(d,r);
	    if(resources_used[t][r]+total > data->getLimit(r,t)) {
	      block_ok=false;  
	      break;
	    }
	  } 
	  if(block_ok) { // use up the resources at t
	    valid_block=true; 
	    int d = block_dest[blk];
	    for(int r = 0; r < r_max; r++){
	      double total = (*blocks)[blk].getRCoef(d,r);
	      resources_used[t][r]+=total;
	    } 
	    // update block time
	    block_times[blk] = t;
	    // now determine the objective value
	    // only if we are within the correct horizon
	    double coef = (*blocks)[blk].getProfit(d);
	    coef /= pow(1+rate,t); 
	    objective += coef;
	    Y[blk][t][d] = 1.0;
	    c1++;
	  }  
	  if(t+1 == horizon && !valid_block) block_times[blk] = -1; //basically didn't complete the block  
	  t++; 
	} 
	bool still_cant = true; 
	for (vector<int>::iterator lit = blocks_waiting.begin(); lit != blocks_waiting.end(); lit++) { 
	  //still_cant = checkPrecedence((*lit),blocks_completed); 
	  //if(!still_cant) { 
	  if(prec_vector[(*lit)] == 0) { 
	    //cout << "Ready to schedule: " << (*lit) << " after: " << blk;
	    //cout << " still_can't is: " << checkPrecedence((*lit),blocks_completed) << endl; 
	    blk=(*lit); 
	    blocks_waiting.erase(lit); 
	    still_cant = false;
	    break;	 
	  } 
	} 
	if(still_cant) blocks_found=false; 
      } 
    } 
  } 
  /*if(blocks_waiting.size()>0){
    cout << "Witing list: " << blocks_waiting.size() << endl;
    exit(0);
  }
  vector<vector<vector<double> > > Yvar(nB);
  for(int b = 0; b < nB; b++){
    Yvar[b].resize(t_max);
    for(int t = 0; t < t_max; t++){
      Yvar[b][t].resize(d_max,0.0);
    }
  }

  int c2 = 0;
  for(int b = 0; b < block_times.size(); b++){
    //if(!bd[b]) continue;
    int blk = block_seq[b];
    if(blk == -1) continue;
    if(block_times[blk] == -1) continue;
    int time = this->block_times[blk];
    int dest = this->block_dest[blk]; // will always be set
    Yvar[blk][time][dest] = 1.0;
    c2++;
    if(Y[blk][time][dest] != 1.0) cout << "Big problem ..." << endl;
  }
  if(c1!=c2){
    cout << "C1 " << c1 << ", C2 " << c2 << endl;
    exit(0);
  }
  for(int r=0; r<r_max; r++){
    for(int t=0; t<t_max; t++){
      double ru = 0.0;
      for(int b=0; b<nB; b++){
	for(int d=0; d<d_max; d++){
	  if(Y[b][t][d] < Yvar[b][t][d])
	    cout << "Yes, this is the problem" << endl;
	  double coef = (*blocks)[b].getRCoef(d,r);
	  ru += coef*Yvar[b][t][d];
	}
      }
      if(ru > data->getLimit(r, t)){
	cout << "Problem at resource " << r << ", time " << t << endl;
	cout << "Available: " << data->getLimit(r, t) << ", used: " << ru << endl;
      }
    }
    }*/
  this->cost = objective;
  return objective;
}


/*
  Set the blocks exactly
*/

double ACO_Solution::setBlocksExact(){
  // set up some parameters
  const int t_max = data->getNPeriod();
  const int r_max = data->getnResources();
  const double rate = data->getDiscountRate();
  std::vector<Block> * blocks=data->getBlock();

  long double objective = 0.0;
  vector<vector<double> > resources_used(t_max);

  for(int t = 0; t < t_max; t++) {
    resources_used[t].resize(r_max,0.0);
  }

  for(int b = 0; b < block_times.size(); b++){
    int t = block_times[b]; 
    if (t == -1) continue;
    int d = block_dest[b];
    for(int r = 0; r < r_max; r++){
      double total = (*blocks)[b].getRCoef(d,r);
      resources_used[t][r]+=total;
      if(resources_used[t][r] > data->getLimit(r,t)) {
	//cout << "Resource problem for block: " << b  << ", used: " 
	//     << resources_used[t][r] << ", available: " << data->getLimit(r,t) << endl;
      }
    } 
    double coef = (*blocks)[b].getProfit(d);
    coef /= pow(1+rate,t); 
    objective += coef;
  } 
  cout << "Inner objective: " << objective << endl;
  this->cost = objective;
  return objective;
}



/*
  Schedule all blocks before considering the waiting list
*/

double ACO_Solution::setBlocksPar(int shift){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int r_max = data->getnResources();
  const double rate = data->getDiscountRate();
  std::vector<Block> * blocks=data->getBlock();

  double objective = 0.0;
  // now set the blocks with precedences
  vector<int> blocks_waiting;

  // set an artificial horizon
  // note, the objective will only include blocks within the original horizon
  int horizon = t_max;
  vector<vector<double> > resources_used(horizon);

  for(int t = 0; t < horizon; t++) {
    resources_used[t].resize(r_max,0.0);
  }
  // determine precedence vector
  vector<int> prec_vector(nB,0);
  for(int b = 0; b < nB; b++){
    prec_vector[b] = (*blocks)[b].getNumPred();    
  }
  // first schedule all blocks
  for(int b = 0; b < block_times.size(); b++){
    int blk = block_seq[b]; // this is needed as blk will change to something on the waiting list
    if(blk == -1) continue;
    // check for successors and negative value
    if(succ[blk].size() == 0 && profit[blk] < 0)  continue;
    if(prec_vector[blk] > 0) blocks_waiting.push_back(blk);
    else{ 
      // update prec_vector for all successors of blk
      for (set<int>::iterator lit = succ[blk].begin(); lit != succ[blk].end(); lit++) { 
	int t_blk = (*lit);
	prec_vector[t_blk]--;
      }

      bool valid_block = false; 
      int t = earliestTime(blk,shift); 
      while(t < horizon && !valid_block){
	bool block_ok = true; 
	// assuming CType = 'L'
	for(int r = 0; r < r_max; r++){
	  double total = 0.0;
	  int d = block_dest[blk];
	  total += (*blocks)[blk].getRCoef(d,r);
	  if(resources_used[t][r]+total > data->getLimit(r,t)) {
	    block_ok=false;  
	    break;
	  }
	} 
	if(block_ok) { // use up the resources at t
	  valid_block=true; 
	  for(int r = 0; r < r_max; r++){
	    double total = 0.0;
	    int d = block_dest[blk];
	    total += (*blocks)[blk].getRCoef(d,r);
	    resources_used[t][r]+=total;
	    // update block time
	    block_times[blk] = t;
	  } 
	  // now determine the objective value
	  // only if we are within the correct horizon
	  int d = block_dest[blk];
	  double coef = (*blocks)[blk].getProfit(d);
	  coef /= pow(1+rate,t); 
	  objective += coef;
	}  
	t++; 
      } 
    } 
  } 

  bool blocks_found = true; 
  int blk = -1;
  for (vector<int>::iterator lit = blocks_waiting.begin(); lit != blocks_waiting.end(); lit++) { 
    if(prec_vector[(*lit)] == 0) { 
      blk=(*lit); 
      blocks_waiting.erase(lit); 
      break;	 
    } 
  } 
  if(blk == -1) cout << "Some big problem ..." << endl;
  while(blocks_found){ 
    // update prec_vector for all successors of blk
    for (set<int>::iterator lit = succ[blk].begin(); lit != succ[blk].end(); lit++) { 
      int t_blk = (*lit);
      prec_vector[t_blk]--;
    }
    bool valid_block = false; 
    int t = earliestTime(blk, shift); 
    while(t < horizon && !valid_block){
      bool block_ok = true; 
      // assuming CType = 'L'
      for(int r = 0; r < r_max; r++){
	double total = 0.0;
	int d = block_dest[blk];
	total += (*blocks)[blk].getRCoef(d,r);
	if(resources_used[t][r]+total > data->getLimit(r,t)) {
	  block_ok=false;  
	  break;
	}
      } 
      if(block_ok) { // use up the resources at t
	valid_block=true; 
	for(int r = 0; r < r_max; r++){
	  double total = 0.0;
	  int d = block_dest[blk];
	  total += (*blocks)[blk].getRCoef(d,r);
	  resources_used[t][r]+=total;
	  // update block time
	  block_times[blk] = t;
	} 
	// now determine the objective value
	// only if we are within the correct horizon
	int d = block_dest[blk];
	double coef = (*blocks)[blk].getProfit(d);
	coef /= pow(1+rate,t); 
	objective += coef;
      }  
      t++; 
    }
    bool cant = true; 
    for (vector<int>::iterator lit = blocks_waiting.begin(); lit != blocks_waiting.end(); lit++) { 
      if(prec_vector[(*lit)] == 0) { 
	blk=(*lit); 
	blocks_waiting.erase(lit); 
	cant = false;
	break;	 
      } 
    } 
    if(cant) blocks_found=false;  
  } 


  // push blocks as early as possible
  // should already be done by the first "time-0" solution

  this->cost = objective;
  return objective;
}


/*
  Set the blocks by time using precedences
  Assumes a sequence of blocks has been defined in block_times
  Assumes each block has been assigned to be sent to a destination
*/

double ACO_Solution::setBlocksWithPrec(const vector<bool> &negBlocks){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  
  // Construct Yvar from block_times and block_dest
  //vector<vector<vector<double> > > Yvar(nB);
  vector<vector<double> > X(nB);
  for(int b = 0; b < nB; b++){
    X[b].resize(t_max,0.0);
    //Yvar[b].resize(t_max);
    //for(int t = 0; t < t_max; t++){
    //  Yvar[b][t].resize(d_max,0.0);
    //}
  }
  for(int b = 0; b < nB; b++){
    int time = block_times[b];
    if( time == -1 ) {
      if(profit[b] < 0){
	time = t_max-1;
      }
      else 
	time = 0;
    }
    int dest = block_dest[b]; // will always be set
    if(time == -1 || dest == -1)
      continue;
    //Yvar[b][time][dest] = 1.0;
    /*if (profit[b] < 0){
      int n_succ = succ[b].size();
      time = max(0,t_max-1-n_succ);
      }*/
    X[b][time] = 1.0;
  }
  // Determines a sequence by sorting by time and then max profit
  Timer timer;
  setSequence(X);
  //cout << "Completed setting sequence, including sorting in " 
  //     << timer.elapsed_time(Timer::REAL) - time_taken << " seconds"
  //     << endl;
  //setSequence(Yvar);
  //double res = setBlocksPar();
  // Do a search for the best negative shift, up to 3 time points
  int shift = 0;
  int limit = 5;
  int best_shift = 0;
  double res = setBlocks(negBlocks, shift);
  while(shift < limit-1){
    shift++;
    double rt = setBlocks(negBlocks, shift);
    if(rt > res){
      res = rt;
      best_shift = shift;
    }
  }
  // finally, use the best_shift
  res = setBlocks(negBlocks, best_shift);
  //cout << "Completed setting blocks in " 
  //     << timer.elapsed_time(Timer::REAL) - time_taken << " seconds"
  //     << endl;
  //cout << "Completed setting blocks ... " << endl;
  return res;
}

/*
  Set the blocks (originally fractional) by time using precedences
  Must generate block_times from the input
  Will always send the blocks to the destination in block_dests
*/

double ACO_Solution::setFracBlocks(const vector<vector<double> > &X_bar, const vector<bool> &negBlocks){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  
  // construct X from X_bar
  vector<vector<double> > X(nB);
  for(int b = 0; b < nB; b++){
    X[b].resize(t_max,0.0);
  }

  for(int b = 0; b < nB; b++){
    double max = -1e30;
    int index = 0;
    for(int t = 0; t < t_max; t++){
      if (X_bar[b][t] > 0.0 && X_bar[b][t] > max){
	max = X_bar[b][t];
	index = t;
      }
    }
    X[b][index] = 1.0;
  }

  // Determines a sequence by sorting by time and then max profit
  setSequence(X);
  int shift = 0;
  return setBlocks(negBlocks, shift);
}

/*
  Given a solution (Yvar or Xvar), try to see if negaitve blocks can be pushed to the next time point
  Sorry, this still needs some thinking. It may only be possible to move a block if it takes up 
    resources from more than one destination. 
*/

double ACO_Solution::postProcess(const vector<vector<vector<double> > > &Y_bar){

}

/*
  swap jobs
 */

ACO_Solution *ACO_Solution::swapBlocks(Random *rnd, const vector<bool> &negBlocks){
  const int nB = block_times.size();
  vector<int> temp_list(block_times); 
  ACO_Solution *new_ant = new ACO_Solution(this->data,this->profit,this->block_dest,this->succ);
  int in_1 = int(double(nB)*rnd->next());
  int in_2 = int(double(nB)*rnd->next());
  
  int j_1 = block_seq[in_1];
  int j_2 = block_seq[in_2];
  temp_list[in_1] = j_2;
  temp_list[in_2] = j_1;
  // now t_list has the new order
  for(int i=0;i<temp_list.size();i++){
    int job = temp_list[i];
    new_ant->block_seq[i] = job;
  }
  new_ant->setBlocksWithPrec(negBlocks);
  return new_ant;
}



/*
  A local search for a solution
  Conduct a systematic 1 swap, always storing an improvement in best
  O(n^2)
 */ 
ACO_Solution* ACO_Solution::localSearch(const vector<bool> &negBlocks){
  ACO_Solution *best = this->copy();
  const int nB = block_times.size();
  int shift = 0;

  for(int i = 0; i < nB; i++){
    for(int j = i+1; j < nB; j++){
      ACO_Solution *ib = best->copy();
      // swap i and j
      int blk = ib->block_seq[i];
      ib->block_seq[i] = ib->block_seq[j];
      ib->block_seq[j] = blk;
      ib->setBlocks(negBlocks, shift);
      if(ib->getObj() > best->getObj()){
	delete best;
	best = ib->copy();
      }
    }
  }
  return best;
}

/*
  Carry out a complex local search
  Not implemented yet, to be done
  DO NOT USE, an empty solution will be returned
*/

ACO_Solution* ACO_Solution::complexLocalSearch(){
  ACO_Solution *t = new ACO_Solution(this->data,this->profit,this->block_dest,this->succ);
  return t;
}

/*
  Beta sampling: take a subsequence of blocks and place them at the
    end of the sequence
*/

ACO_Solution *ACO_Solution::betaSampling(Random *rnd, const vector<bool> &negBlocks){
  const int nB = block_times.size();
  vector<int> temp_list;
  vector<int> t_list;
  ACO_Solution *new_ant = new ACO_Solution(this->data,this->profit,this->block_dest,this->succ);

  int s = double(nB)*rnd->next();
  int e = double(nB-s)*rnd->next();
  
  for(int i = s; i < e+s ; i++) temp_list.push_back(block_seq[i]);
  for(int i = 0; i < s; i++) { // fill up all the jobs until s
    int job = this->block_seq[i];
    new_ant->block_seq[i] = job;
  }
  for(int i = e+s ; i < nB; i++){ // move the jobs from the end to s in the new schedule
    int job = this->block_seq[i];
    new_ant->block_seq[i] = job;
  }
  double beta_s = 1.0-(5.0/(double)temp_list.size());
  int try1 =0;
  while(temp_list.size()>0){
    for(int j=0;j<temp_list.size();j++){
      double p = rnd->next();
      if(p >= beta_s) try1 ++;
      else {
	t_list.push_back(temp_list[j]);
	try1 = 0;
      }
    }
    
    if(try1 >= 100) {
      for(int k=0;k<temp_list.size();k++) t_list.push_back(temp_list[k]);
    }
    for(int j=0;j<t_list.size();j++) {
      int job = t_list[j];
      for (vector<int>::iterator lit = temp_list.begin(); lit != temp_list.end(); lit++) {
	if((*lit)==job) {
	  temp_list.erase(lit);
	  break;
	}
      }
    }
  }
  // now t_list has the new order
  for(int i=s;i<e+s;i++){
    int blk = t_list[i-s];
    new_ant->block_seq[i] = blk;
  }
  int shift = 0;
  new_ant->setBlocks(negBlocks, shift);
  return new_ant;
}

 
/* 
   Function to select times for blocks
   Also returns a destimation for a block
*/ 
pair<int,int> ACO_Solution::bestSelect(Pheromones *p, int variable){ 
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();

  // Determine the time and destination
  double value_best = 0.0; 
  int time = -1;
  int dest = -1;
  for ( int l = 0 ; l < t_max+1 ; l++ ) { 
    double heuristic = p->get_ci(variable, l); 
    if ( heuristic > value_best ) { 
      value_best = heuristic; 
      time = l; 
    } 
  } 
  if(value_best == 0.0 || time == t_max) {
    return make_pair(-1,-1);
  }  

  value_best = 0.0;
  for ( int l = 0 ; l < d_max ; l++ ) { 
    double heuristic = p->get_ci_y(variable, l); 
    if ( heuristic > value_best ) { 
      value_best = heuristic; 
      dest = l; 
    } 
  }  

  if(value_best == 0.0) {
    return make_pair(-1,-1);
  } 

  return make_pair(time,dest); 	 
} 
 
pair<int,int> ACO_Solution::probSelect(Pheromones *p, int variable, Random *rnd){ 
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();

  double sumProbs=0.0; 
  double prob=0.0; 
  vector<double> selectionProb_time(t_max,0.0); 
  vector<double> selectionProb_dest(d_max,0.0); 
  
  for(int j = 0; j < t_max; j++){ 
    selectionProb_time[j] = p->get_ci(variable, j); 
    sumProbs += selectionProb_time[j];
  }	 
  // no block can be selected
  if(sumProbs == 0.0){
    return make_pair(-1,-1);
  }

  double r_num= rnd->next();   
  double r = sumProbs*r_num;	 
  int time = 0; 
  prob = selectionProb_time[time]; 
  while( prob < r ){ 
    time++; 
    if(time == t_max){
      break;
      cout << "Time problem ..." << endl;
      exit(0);
    }
    prob += selectionProb_time[time]; 
  }	 
  // chosen not to mine the block
  if(time == t_max){
    return make_pair(-1,-1);
  }
  //cout << "Selected time: " << time << endl;
  sumProbs = 0.0;
  for(int j = 0; j < d_max; j++){ 
    selectionProb_dest[j] = p->get_ci_y(variable, j); 
    sumProbs += selectionProb_dest[j]; 
  }	 

  // no block can be selected
  if(sumProbs == 0.0){
    return make_pair(-1,-1);
  }

  r_num= rnd->next();   
  r = sumProbs*r_num;	 
  int dest = 0; 
  prob = selectionProb_dest[dest]; 
  while( prob < r ){ 
    dest++; 
    if(dest == d_max){
      break;
      cout << "Destination problem ..." << endl;
      exit(0);
    }
    prob += selectionProb_dest[dest]; 
  }	 
  return make_pair(time, dest);	 
} 
 
 
/* 
   Select a block in the usual ACO way
     Uses greedy selection or heuristic selection
   Earliest is specifically for selecting a time
*/ 
pair<int,int> ACO_Solution::selectBlock(int variable, Pheromones *p, Random *r, double q_0){ 
  pair<int,int> block_dest;
  double q=r->next(); 
  if(q<q_0) 
    block_dest = probSelect(p, variable, r); 
  else 
    block_dest = bestSelect(p, variable); 
  
  // set the values
  this->block_times[variable] = block_dest.first;
  //this->block_dest[variable] = block_dest.second;

  return block_dest; 
} 

/*
  Get the times assigned to each block
*/

const vector<int> & ACO_Solution::getBlockTimes(){
  return this->block_times;
}

/*
  Get the destinations assigned to each block
*/

const vector<int> & ACO_Solution::getBlockDests(){
  return this->block_dest;
}

/*
  Get the time assigned to a block
  Note: could be -1, i.e. the block is not assigned
*/
int ACO_Solution::getBlockTime(int block){
  return block_times[block];
}

/*
  Get the destination assigned to a block
  Note: could be -1, i.e. the block is not assigned
*/

int ACO_Solution::getBlockDest(int block){
  return block_dest[block];
}

/*
  return the objective value
*/

double ACO_Solution::getObj(){
  return this->cost;
}

/*
  return the cumulative violation
*/

double ACO_Solution::getVio(){
  return this->violations;
}

/*
  Return a Y var solution
*/

vector<vector<vector<double> > > ACO_Solution::mappedYvar(){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const int r_max = data->getnResources();
  std::vector<Block> * blocks=data->getBlock();

  vector<vector<vector<double> > > Yvar(nB);
  for(int b = 0; b < nB; b++){
    Yvar[b].resize(t_max);
    for(int t = 0; t < t_max; t++){
      Yvar[b][t].resize(d_max,0.0);
    }
  }

  for(int b = 0; b < block_times.size(); b++){
    int blk = block_seq[b];
    if(blk == -1) continue;
    if(block_times[blk] == -1) continue;
    int time = this->block_times[blk];
    int dest = this->block_dest[blk]; // will always be set
    Yvar[blk][time][dest] = 1.0;
  }
  for(int r=0; r<r_max; r++){
    for(int t=0; t<t_max; t++){
      double ru = 0.0;
      for(int b=0; b<nB; b++){
	for(int d=0; d<d_max; d++){
	  double coef = (*blocks)[b].getRCoef(d,r);
	  ru += coef*Yvar[b][t][d];
	}
      }
      if(ru > data->getLimit(r, t)){
	cout << "Problem at resource " << r << ", time " << t << endl;
	cout << "Available: " << data->getLimit(r, t) << ", used: " << ru << endl;
      }
    }
  }
  return Yvar;
}

/*
  Update the profit vector
*/

void ACO_Solution::updateProfit(const vector<double> &uprofit){
  this->profit.clear();
  this->profit = uprofit;
}

/*
  New set of functions for second ACO scheme
  The new scheme does not worry about resources and chooses
    a block at a time point
  After computing the NPV it computes the violations:
    I.e. the amount of resource violated at every time point.
  Assumes a sequence is defined, can be -1, i.e. a block is not used
 */

double ACO_Solution::computeNPV(){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const int r_max = data->getnResources();
  const double rate = data->getDiscountRate();
  std::vector<Block> * blocks=data->getBlock();

  // resources used per time point
  vector<vector<int> > ru(t_max);

  for(int t = 0; t < t_max; t++) {
    ru[t].resize(r_max,0.0);
  } 
  double objective = 0;

  for (int b = 0; b < nB; b++){
    int t = block_times[b];
    if (t == -1) // unscheduled block 
      continue;
    int dest = block_dest[b];
    double coef = (*blocks)[b].getProfit(dest);
    coef /= pow(1+rate,t); 
    objective += coef;
    // compute the resources used
    for(int r = 0; r < r_max; r++){
      double total = 0.0;
      for(int d=0; d<d_max; d++){ 
	double coef = (*blocks)[b].getRCoef(d,r);
	total += coef;
      }
      ru[t][r] += total;
    }
  }
  // Set the solution's cost
  cost = objective;

  // Compute violations at every time point
  violations = 0.0;
  for(int t = 0 ; t < t_max; t++){
    for(int r = 0; r < r_max; r++){
      if(ru[t][r] > data->getLimit(r,t)){
	violations += ru[t][r] - data->getLimit(r,t);
      }
    }
  }
  return objective;
}

/*
  Compare two solutions using violations as the first measure
*/

double ACO_Solution::compareSol(ACO_Solution *a){
  return (this->violations < a->violations ||
	  (this->violations == a->violations && this->cost > a->cost));
}

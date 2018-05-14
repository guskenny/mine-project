/*************************************************************************** 
                            ACO_solution.h 
                         ------------------- 
    last modified   : 13/5/2016 
    copyright       : (C) 2016 by Dhananjay Thiruvady 
    libraries		: . 
    description		: a solution to the problem as viewed by ACO 
 ***************************************************************************/ 
 
/*************************************************************************** 
 *                                                                         * 
 *   This program is free software; you can redistribute it and/or modify  * 
 *   it under the terms of the GNU General Public License as published by  * 
 *   the Free Software Foundation; either version 2 of the License, or     * 
 *   (at your option) any later version.                                   * 
 *                                                                         * 
 ***************************************************************************/ 

#ifndef ACO_Solution_H
#define ACO_Solution_H

#include <iostream>
#include <new>
#include <vector>
#include <set>
#include <utility>

using namespace std;

#include "Random.h"
#include "pheromones.h"

class ACO_Solution{
 private:
  double cost;
  double violations;

  // vector of times for each block
  vector<int> block_times;
  // determine a sequence of the blocks
  vector<int> block_seq;
  // vector of destinations for each block
  // Currently (6/7/16), blocks can only be sent to one destination
  // This should change to proportions, but needs some thinking
  vector<int> block_dest;

  vector<double> profit;
  Daten *data;

  // list of successors for each block
  vector<set<int> > succ;

  // Some local functions
  void setProfit();
  bool checkPrecedence(int blk, const vector<bool> &blocks_completed);
  ACO_Solution *swapBlocks(Random *rnd, const vector<bool> &negBlocks);
  pair<int,int> bestSelect(Pheromones *p, int variable);
  pair<int,int> probSelect(Pheromones *p, int variable, Random *rnd);
  int earliestTime(int b, int shift);
  void sortSeq(vector<vector<int> > &timeBlocks);
  void sortProbSeq(vector<vector<int> > &timeBlocks);
  
 public:
  //ACO_Solution(Daten *data);
  ACO_Solution(Daten *data, const vector<double> &profit, 
		 const vector<int> &bd, const vector<set<int> > &succ);

  ACO_Solution(Daten *data, const vector<double> &profit, 
	       const vector<int> &bd, const vector<set<int> > &succ,
	       const vector<bool> &negBlocks,
	       const vector<int> &X,
	       const vector<vector<double> > &Y,
	       bool exact
	       ); 

  ACO_Solution *copy();
  pair<int,int> selectBlock(int variable, Pheromones *p, Random *r, double q_0);

  ACO_Solution* localSearch(const vector<bool> &negBlocks);
  ACO_Solution* complexLocalSearch();
  ACO_Solution* betaSampling(Random *rnd, const vector<bool> &negBlocks);


  double setBlocks(const vector<bool> &negBlocks, int shift);
  double setBlocksExact();
  double setBlocksPar(int shift);
  double setBlocksWithPrec(const vector<bool> &negBlocks);
  double setFracBlocks(const vector<vector<double> > &X_bar, const vector<bool> &negBlocks);
  void setTimes(const vector<vector<vector<double> > > &Yvar);
  void setTimes(const vector<int> &times);
  void setSequence(const vector<vector<vector<double> > > &Yvar);
  void setSequence(const vector<vector<double> > &X);
  void setBackPropProfit();
  double postProcess(const vector<vector<vector<double> > > &Y_bar);


  // retrieve some information
  const vector<int> & getBlockTimes();
  const vector<int> & getBlockDests();
  int getBlockTime(int block);
  int getBlockDest(int block);
  double getObj();
  double getVio();
  vector<vector<vector<double> > > mappedYvar();

  void updateProfit(const vector<double> &uprofit);

  // New set of functions with second ACO scheme
  double computeNPV();
  double compareSol(ACO_Solution *a);
  ~ACO_Solution(){};
};

#endif

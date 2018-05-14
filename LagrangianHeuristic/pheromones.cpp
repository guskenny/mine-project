/***************************************************************************
                            pheromones.cpp
                         -------------------
    last modified   : 07/07/2016
    copyright       : (C) 2008 by Dhananjay Thiruvady
    libraries		: .
    description		: The pheromones class
                          A pheromone T_ab defines the probability of 
                          selecting block a at time point b 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "pheromones.h"

Pheromones::Pheromones(Daten *data, const vector<double> & profit){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const double rate = data->getDiscountRate();
  //std::vector<Block> * blocks=data->getBlock();

  p.resize(nB);
  ci.resize(nB);
  p_y.resize(nB);
  ci_y.resize(nB);
  for(int i = 0; i < nB; i++){
    p[i].resize(t_max+1,1.0/double(t_max));
    ci[i].resize(t_max+1,1.0/double(t_max));
    p_y[i].resize(d_max,1.0/double(d_max));
    ci_y[i].resize(d_max,1.0/double(d_max));
  }		
  tau_min = 0.001;
  tau_max = 0.999;
  alp = 1.0; // probably not optimal
  bet = 1.0; // probably not optimal


  // determine the profit matrix for the heuristic information
  // First determine the smallest value
  double minProf = 1e99;
  double maxProf = -1e99;
  for(int b = 0; b < nB; b++){
    if(profit[b] < minProf) minProf = profit[b];
    if(profit[b] > maxProf) maxProf = profit[b];
  }
  double divFactor = 1.0;
  while(maxProf > 1.0) {
    maxProf /= 10.0;
    divFactor *= 10.0;
  }
  heuristicProfitMatrix.resize(nB);
  for(int b = 0; b < nB; b++){
    heuristicProfitMatrix[b].resize(t_max+1,0.0);
    for(int t = 0; t < t_max; t++){
      double factor = (profit[b]+abs(minProf))/divFactor;
      if(profit[b] > 0.0){
	factor /= pow(t+1,2);
      }
      else {
	factor /= 10.0;
	factor *= pow(t+1,2);
      }
      //cout << "b: " << b << ", t:" << t << ", factor: " << factor << ", profit: " << profit[b] << endl;
      heuristicProfitMatrix[b][t] = factor;
    }
  }
}

/*
  Reset the pheromone values to their original state
  Could have done this by deleting objects ... may implement it later
*/

void Pheromones::reset(Daten *data){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();

  for(int i = 0; i < nB; i++){
    for(int j = 0; j < t_max; j++){
      p[i][j] = 1.0/double(t_max);
      ci[i][j] = 1.0/double(t_max);
    }
    for (int j = 0; j < d_max; j++){
      p_y[i][j] = 1.0/double(d_max);
      ci_y[i][j] = 1.0/double(d_max);
    }
    p[i][t_max] = 1.0/double(t_max);
    ci[i][t_max] = 1.0/double(t_max);
  }		
}

// accessor methods
double Pheromones::get_p(int j1, int j2){
  return p[j1][j2];
}
double Pheromones::get_ci(int j1, int j2){
  return ci[j1][j2];
}
double Pheromones::get_p_y(int j1, int j2){
  return p_y[j1][j2];
}
double Pheromones::get_ci_y(int j1, int j2){
  return ci_y[j1][j2];
}

//mutator methods
void Pheromones::set_p(int j1, int j2, double val){
  p[j1][j2]=val;
}
void Pheromones::set_ci(int j1, int j2, double val){
  ci[j1][j2]=val;
}

/*
  Local pheromone update for an ACS style algorithm 
*/

void Pheromones::localPheromoneUpdate(int block, int time, int dest, double lrate){ 
  int t = time;
  if(time == -1) t = p[block].size()-1;
  this->p[block][t] = (1.0 - lrate) * this->p[block][t] + tau_min;
  // determine the heuristic component based on a block's profit
  //cout << t << ":" << d << ", size: " << heuristicProfitMatrix[block].size() <<  endl;
  double heuristic = heuristicProfitMatrix[block][t];  
  // the heuristic needs some thinking
  this->ci[block][t] =  pow(this->p[block][t],alp) * pow(heuristic,bet); 
  if(dest != -1){ // there is no destination if we don't select a block
    this->p_y[block][dest] = (1.0 - lrate) * this->p[block][dest] + tau_min;
    // the heuristic needs some thinking
    this->ci_y[block][dest] =  pow(this->p_y[block][dest],alp) * pow(heuristic,bet); 
  }
}

/*
  Global pheromone update for ACS
*/

void Pheromones::globalPheromoneUpdate(const vector<int> &blockTimes, const vector<int> &blockDests, double lrate, double obj){
  for ( int i = 0 ; i < p.size(); i++ ) {
    for(int j = 0; j < p[i].size(); j++){
      p[i][j] = (1.0-lrate)*p[i][j];
      if(p[i][j] < tau_min) 
	p[i][j] = tau_min;
      // the heuristic needs some thinking
      double heuristic = heuristicProfitMatrix[i][j];  
      ci[i][j] = pow(p[i][j],alp) * pow(heuristic,bet); 
    }    
  }

  double delta_tau = abs(obj);
  while(delta_tau > 10.0) delta_tau /= 10.0;
  //cout << "DT: " << delta_tau << endl;
  const int nB = blockTimes.size(); // should be the number of blocks

  for ( int i = 0 ; i < nB; i++ ) {
    int j = blockTimes[i];
    int k = blockDests[i];
    bool skip = false;
    if(j == -1){
      j = heuristicProfitMatrix[i].size()-1;
      k = 0;
      skip = true;
    }
    double heuristic = heuristicProfitMatrix[i][j];  
    //delta_tau = heuristic;
    //p[i][j] = (1.0-lrate)*p[i][j] + delta_tau;
    p[i][j] += delta_tau;
    if(p[i][j] > tau_max) 
      p[i][j] = tau_max;
    // the heuristic needs some thinking
    ci[i][j] = pow(p[i][j],alp)* pow(heuristic,bet); 
    //cout << i << ":" << j << ", p: " << p[i][j] << ", ci: " << ci[i][j] <<  endl;
    if(!skip){
      p_y[i][k] = (1.0-lrate)*p_y[i][k] + delta_tau;
      if(p_y[i][k] > tau_max) 
	p_y[i][k] = tau_max;
      if(p_y[i][k] < tau_min) 
	p_y[i][k] = tau_min;
      // the heuristic needs some thinking
      ci_y[i][k] = pow(p_y[i][k],alp);//* pow(heuristic,bet);     
    }
  }
}

/*
  A simple function to display the pheromones
  Mainly for validation
*/

void Pheromones::display(){
  cout << "Pheromones: " << endl;
  for(int i = 0; i < p.size(); i++){
    for(int j = 0; j < p[i].size(); j++){
      cout << "(" << i << "," << j << ")" << p[i][j] << ", ";
    }
    cout  << endl;
  }
  cout << "Choice information: " << endl;
  for(int i = 0; i < p.size(); i++){
    for(int j = 0; j < p[i].size(); j++){
      cout << "(" << i << "," << j << ")" << p[i][j] << ", ";
    }
    cout << endl;
  }
}

/***************************************************************************
                            pheromones.h
                         -------------------
    last modified   : 1/5/2008
    copyright       : (C) 2008 by Dhananjay Thiruvady
    libraries		: .
    description		: The header for the Pheromone class
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 
#ifndef Pheromones_H
#define Pheromones_H

#include <iostream>
#include <new>
#include <vector>
#include <cmath>

using namespace std;

#include "../include/daten.h"

class Pheromones{
 private:
  vector<vector<double> >  p;
  vector<vector<double> >  ci;
  vector<vector<double> >  p_y;
  vector<vector<double> >  ci_y;
  
  // Some constants specific to pheromones
  // At a future point, may want to enter these as parameters
  double tau_min;
  double tau_max;
  double bet; // probably not optimal
  double alp; // probably not optimal
  
  // keep a profit matrix used for heuristic information
  vector<vector<double > > heuristicProfitMatrix;

  // Some local functions

  // mutator functions
  void set_p(int j1, int j2, double val);
  void set_ci(int j1, int j2, double val);

 public:
  
  Pheromones(Daten *data, const vector<double> &profit);
  void reset(Daten *data);
  // Pheromone update functions
  void localPheromoneUpdate(int block, int time, int dest, double lrate);
  void globalPheromoneUpdate(const vector<int> &blockTimes, const vector<int> &blockDests, double lrate, double obj);
  void display();
  
  // accessor functions
  double get_p(int j1, int j2);
  double get_ci(int j1, int j2);
  double get_p_y(int j1, int j2);  
  double get_ci_y(int j1, int j2);  
  ~Pheromones(){};

};

#endif

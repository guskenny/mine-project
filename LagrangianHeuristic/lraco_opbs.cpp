/***************************************************************************
                            Lagrangian relaxation for the OPBS problem
                         -------------------
    last modified   : 13/5/2016
    copyright       : (C) 2016 by Dhananjay Thiruvady
    description	    : main file to get the algorithms going
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
#include <string>
#include <cmath>
#include <cstring>
#include <iomanip>


using namespace std;

#include "Random.h"
#include "Timer.h"
#include "../include/daten.h"
#include "solver.h"
#include "solver_functions.h"

/*
  the following variables could be modified by input from the command line
*/


// count the number of labeling steps, set maximum to 1 trillion
long int steps = 100000000;
bool steps_given=false;

// variables to test parameters
bool time_limit_given = false;
bool iter_limit_given = false;
bool input_file_given = false;
string input_file;
double time_limit=900.0;
bool use_prec=false;

// display information through the run or not?
bool disp = true;
// show the state of the solver
bool show = true;

// buffer for deadlines, only for Vanhouke 2010: start with the tightest deadline
int buffer = 5;

// only run one of the following: 0 = MIP, 1 = LR, 2 = NF 
int alg_type = 0; // using the lagrangian relaxation scheme

// which MIP model, simple for original
bool simple = true;

// which network flow model? Primal by default 
int nf_type = 0;

// provide a bound which will improve the subgradient optimisation
double bound = -1.0;

// ACO parameters
int nants = 10;
double lrate = 0.01;
double q_0 = 0.5;
int max_iter = 1000;

/* BeamACO parameters */
int beam_width = 10;
int mu = 2;

// use ACO within LR or not, not by default
bool use_aco = false;

// Tolerance for difference amount
int diff_amt = 1e2;

// read the parameters provided as input
void read_parameters(int argc, char **argv){
  int iarg=1;
  while (iarg < argc)  {
    if (strcmp(argv[iarg],"-f")==0) {
      input_file = argv[++iarg];
      input_file_given = true;
    }
    else if (strcmp(argv[iarg],"-t")==0) { 
      time_limit=atof(argv[++iarg]);
      time_limit_given = true;
    }
    else if (strcmp(argv[iarg],"-buffer")==0) { 
      buffer=atoi(argv[++iarg]);
    }
    else if (strcmp(argv[iarg],"-maxiter")==0) {
      max_iter=atoi(argv[++iarg]);
      iter_limit_given = true;
    }
    else if (strcmp(argv[iarg],"-nants")==0) {
      nants = atoi(argv[++iarg]);
    }
    else if (strcmp(argv[iarg],"-w")==0) {
      beam_width = atoi(argv[++iarg]);
    }
    else if (strcmp(argv[iarg],"-lrate")==0) {
      lrate = atof(argv[++iarg]);
    }
    else if (strcmp(argv[iarg],"-q_0")==0) {
      q_0 = atof(argv[++iarg]);
    }
    else if (strcmp(argv[iarg],"-steps")==0) {
       steps = atoi(argv[++iarg]);
	steps_given = true;
    }
    else if (strcmp(argv[iarg],"-mu")==0){  
      mu=atoi(argv[++iarg]);
    }
    else if (strcmp(argv[iarg],"-disp")==0){  
      disp=atoi(argv[++iarg]);
    }
    else if (strcmp(argv[iarg],"-show")==0){  
      show=atoi(argv[++iarg]);
    }
    else if (strcmp(argv[iarg],"-rho")==0){  
      lrate=atof(argv[++iarg]);
    }
    else if (strcmp(argv[iarg],"-alg")==0){  
      alg_type=atoi(argv[++iarg]);
    }
    else if (strcmp(argv[iarg],"-simple")==0){  
      simple=atoi(argv[++iarg]);
    }
    else if (strcmp(argv[iarg],"-nf")==0){  
      nf_type=atoi(argv[++iarg]);
    }
    else if (strcmp(argv[iarg],"-bound")==0){  
      bound=atof(argv[++iarg]);
    }
    else if (strcmp(argv[iarg],"-useaco")==0){  
      use_aco=atoi(argv[++iarg]);
    }
    else if (strcmp(argv[iarg],"-diffamt")==0){  
      diff_amt=atoi(argv[++iarg]);
    }
    iarg++;
  }
  if (input_file_given == false) {
    printf("No input file?");
    printf("\n");
    exit(1);
  }
}

/* controller method */

int main(int argc, char *argv[]){
  ifstream in,in2;
  time_t t1;
  long ltime;
  int stime;
  
  ltime=time(NULL);
  stime=(unsigned) ltime/2;
  srand(ltime);
  
  (void) time(&t1);
  if ( argc < 2 ) {
    cout << "\nError with input\n"; 
    exit(1);
  }
  else {
    read_parameters(argc,argv);
    in.open(input_file.c_str());
    if(!in) { // file couldn't be opened
      cout << "Error: file does not exist" << endl;
      exit(1);
    }
    Daten *data = new Daten(input_file.c_str(),'p');
    const char pT = data->getProbType();
    const int nB = data->getNBlock();
    const int t_max = data->getNPeriod();
    const int d_max = data->getnDestination();
    const int r_max = data->getnResources();
    const double rate = data->getDiscountRate();
    std::vector<Block> * blocks=data->getBlock();
    cout<<setprecision(12);
    // info on screen
    {
      std::string datafile(argv[2]);
      std::string strng;
      std::string::size_type foundPos, foundEnd;
      foundPos = datafile.find_last_of("/");
      if(foundPos==std::string::npos) strng="";
      else strng = datafile.substr(foundPos+1);
      cout<<"Data: " << strng << endl;
      cout<<"Type: ";
      if(pT=='p') cout<< "PCPSP"<<endl;
      else if (pT=='c') cout<< "CPIP"<<endl;
      else cout<< "UPIP"<<endl;
      cout<<"NBlocks: "<< nB<<endl;
      if(pT!='u'){
	cout<<"NPeriods: "<< t_max<<endl;
	cout<<"NResource_Side_Constraints: "<< r_max<<endl;
	cout<<"Discount_Rate: "<< rate<<endl;	 
	if(pT!='c'){
	  cout<<"NDestination: "<< d_max<<endl;
	}
      }
    }
    // Create a MaxClosure object
    SolverFunctions *funcs = new SolverFunctions(*data, diff_amt);
    Solver *alg_obj = new Solver(*data, *funcs);
    
    // simple is only valid for the mip
    // nf is only valid for the nf algorithm
    alg_obj->run_algorithm(disp, time_limit, show, alg_type, 
			   simple, nf_type, bound, nants, lrate, q_0, use_aco);
    delete funcs;
    delete alg_obj;
    delete data;
  }
  in.close();
  return 0;
}

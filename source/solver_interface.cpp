/******************************************************************************** 
                            Solver Interface
                         -------------------------- 
    last modified   : 25/08/2016 
    copyright       : (C) 2016 by Dhananjay Thiruvady 
    libraries: . 
    description: Implementation of the functions in solver interface
********************************************************************************/ 
 
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
#include <cmath>

using namespace std;

#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp> // define boost::tie
#include "MaxClosure_Base.h"
#include "MaxClosure_BoostMaxFlow_BK.h"
#include "MaxClosure_BoostMaxFlow_EK.h"
#include "MaxClosure_BoostMaxFlow_PR.h"
#include "MaxClosure_NetworkFlow.h"
#include "graph.h"
#include "../QOL/CpuTimer.h"
#include "MaxClosureFactory.h"
#include "Particle.h"// LaPSO & Volume algorithms
#include "bz.h"
#include "solver_interface.h"
#include "Preprocess.h"

SolverInterface::SolverInterface(Daten &data){
  this->data = new CumulativeModel(data);
  deleteData=true;
  mip_sol = NULL;
}

SolverInterface::SolverInterface(SolverInterface &other)
{
  if(other.deleteData)
    data = new CumulativeModel(*other.data);
  else
    data = other.data;
  deleteData = other.deleteData;
  mip_sol = NULL;
}

// Generate a Y variable like input for the algorithms
void SolverInterface::convertToY(const BranchNode_info &bni, vector<vector<vector<double> > > &Y){

}

// Convert the current solution to an integer solution
Sol_Int *SolverInterface::convertToSol_Int(const vector<int> &X, 
					   const vector<vector<double> > &Y, 
					   const double &obj
					   ){
  const int t_max = data->getNPeriod();
  Sol_Int *sol = new Sol_Int();
  sol->x = X;
  sol->y = Y;
  sol->nT = t_max;
  sol->obj = obj;
  return sol;
}

/*
  Convert data from an ACO solution into vectors needed by problemData

*/

void SolverInterface::collecDataFromACO(ACO_Solution *a, vector<int> &X, vector<vector<double> > &Y){
  const int t_max = data->getNPeriod();
  // Do the conversions
  const int d_max = data->getnDestination();  
  vector<int> times = a->getBlockTimes();
  vector<int> dests = a->getBlockDests(); // this is limiting in the case of the MIP
  
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
}

void SolverInterface::storeCumulativeSoln(Sol_Real &sol,const vector<double> &y)
{
  sol.init(data->getNBlock(),data->getNPeriod(),data->getnDestination());
  for(int b=0;b<data->getNBlock();++b)
    for(int t=0;t<data->getNPeriod();++t){
      sol.x[b][t]=0;
      for(int d=0;d<data->getnDestination();++d)
	sol.x[b][t]+= sol.y[b][t][d]
	  = data->blockProcessed(y,b,d,t);
    }
}
void SolverInterface::storeCumulativeSoln(Sol_Int &sol,const vector<double> &y)
{
  sol.init(data->getNBlock(),data->getnDestination());
  sol.nT = 0;
  sol.obj=0;
  for(int b=0;b<data->getNBlock();++b)
    for(int t=0;t<data->getNPeriod();++t){
      if(data->blockMined(y,b,t)){
	sol.x[b]=t;
	if(t>=sol.nT) sol.nT=t+1;
	for(int d=0;d<data->getnDestination();++d){
	  sol.y[b][d]=data->blockProcessed(y,b,d,t);
	  sol.obj+=data->getProfit(b,d,t);
	}
	break;
      }
    }
}

void SolverInterface::storeCumulativeSoln(Sol_Int &sol,const vector<int> &x)
{   // convert x to double vector - crude but effective
  std::vector<double> y(x.size());
  for(size_t i=0;i<x.size();++i) y[i] = x[i];
  storeCumulativeSoln(sol,y);
}


/* 
   Call one of the solvers (e.g. Lagrangian relaxation) using the Y variable
   Returns success (0) or failure (1), or maybe something else?
   alg_type defines which algorithm is going to be called.
     (1) ACO
     (2) BZ // LR (no ACO)
     (3) LR (with ACO)
     (4) MIP-ACO
     (5) BZ-MIP-ACO? -> need to implement
     (6) Volume
     (7) LaPSO

   Note some of these will call others as part of what they do, e.g. Lagrangian relaxation 
      will make use of ACO
*/

int SolverInterface::callSolver(MineProblem &problem, int diff_amt, int iterations, int alg_type, int time, int shift){
  //cout << "Running heuristics ... " << endl;

  bool disp = true;
  bool show = true;

  Preprocess preprocess(*problem.data,problem.get_BranchNode());
  std::cout << "After preprocessing, branch has " << problem.get_BranchNode().fixedCnt() << " fixed\n";
  SolverFunctions *funcs = new SolverFunctions(*data, diff_amt);
  Solver *alg_obj = new Solver(*data, *funcs);
  // Construct simple solution from problem info 
  // Note, it won't be exact, since the heuristic will be used to re-order the blocks
  //vector<int> X = problem.current->sol_int.x;
  //vector<vector<double> > Y = problem.current->sol_int.y;
  vector<int> X;
  vector<vector<double> > Y;
  int nT;
  double obj;
  problem.get_sol(nT, obj, X, Y);
  obj = -1e30;
  const int t_max = data->getNPeriod();
  const int nB = data->getNBlock();
  // Construct a solution from problem or X and Y 
  ACO_Solution* sol = new ACO_Solution(data,funcs->profit,
				       funcs->block_dest,funcs->succ,
				       funcs->negBlocks,
				       X, Y, false); 
  char maxClosureMethod='B'; double timeLimit=3600; // REPLACE THESE WITH ACTUAL PARAMETERS
  MaxClosureFactory mcfactory(maxClosureMethod);
  if(alg_type == 1){ // ACO, 2 needs to change to the BZ algorithm
    Pheromones *p = new Pheromones(data, funcs->profit);
    // some obvious parameters have been hard coded
    ACO_Solution* a = alg_obj->ACO_run(p, sol, iterations, 20, 0.2, 0.1, disp, 1000, &problem);
    collecDataFromACO(a, X, Y);
    obj = a->getObj();
    cout << "\tCompleted ACO with objective: " << obj << endl;
    //cout << "best MIP solution: " << mip_obj << endl;
    delete a;  
    delete p;
  }else if( alg_type == 2 ){ // BZ algorithm
    BZ bz(*data,mcfactory,timeLimit,&problem.get_BranchNode()); 
    bz.solve(&problem);		// best integer solution updated during run
    Sol_Real sol;
    problem.get_sol(sol);
    std::cout << "BZ method terminated with UB="<< sol.obj<<"using " << sol.nT << " periods\n";
  }
  else if(alg_type == 3){ // LR with ACO, change the last arguement to 0 to not use ACO.
    bool time_limit = iterations;
    double bound = 0; // really does nothing, was used to get nicer multipliers 
    ACO_Solution* a= alg_obj->Lagrangian_heuristic(false, time_limit, false,  bound, iterations, 20, 0.2, 0.1, 1, &problem);
    collecDataFromACO(a, X, Y);
    obj = a->getObj();
    cout << "\tCompleted LR-ACO with objective: " << obj << endl;
    delete a;
  }
  else if(alg_type == 4){
    // do the VLNS
    int time_limit = iterations; // an approximation for now
    // X and Y will contain the solution information when done
    obj = alg_obj->LargeNeighbourhoodSearch( disp, time_limit, show, sol, X, Y, &problem);
    cout << "\tCompleted VLNS with objective: " << obj << endl;
  }
  else if(alg_type == 5){
    cout << "\tStarted large neighbourhood search, fixed ... " << endl;
    // do the VLNS
    int time_limit = iterations; // an approximation for now
    // fix some dummy node info, i.e., every task may only be done at the first time_point
    for(int b = 0; b < nB; b++){
      std::array<int, 2> tim = {0,1};
      problem.get_BranchNode().time.push_back(tim);
      // the statements below are not needed, just testing
      problem.get_BranchNode().time[b][0]=0;
      problem.get_BranchNode().time[b][1]=1;
    }
    cout << "\tFinished determining BNI ... " << endl;
    // X and Y will contain the solution information when done
    obj = alg_obj->LNS_fixed( disp, time_limit, show, &problem, X, Y);
    cout << "\tCompleted VLNS with objective: " << obj << endl;
  }else if(alg_type==6 || alg_type==7 ) { // Volume algorithm or LaPSO
    std::vector<const char *> arg;
    if(alg_type ==6)arg.push_back("-V"); else arg.push_back("-L");
    arg.push_back("--maxIter");
    std::string maxIter=boost::str(boost::format("%d")%iterations);
    arg.push_back(maxIter.c_str());
    arg.push_back("--maxWall");
    std::string maxWall=boost::str(boost::format("%.2f")%timeLimit);
    arg.push_back(maxWall.c_str());
    LaPSO::Problem *lapso = runParticleMain(*data,arg.size(),&arg[0]);
    LaPSO::Particle &best=lapso->best;
    // Note: costs are -ve, pertub = fractiona, x = interger solution
    Sol_Real solUB(data->getNBlock(),data->getNPeriod(),data->getnDestination());
    Sol_Int solLB(data->getNBlock(),data->getnDestination());
    data->storeSoln(solUB,best.perturb);
    data->storeSoln(solLB,best.x); // sets .obj correctly
    std::cout << "Vol/Particle found " << solLB.obj << " " << solUB.obj<<std::endl;
    problem.update_sol(solUB);
    problem.update_sol(solLB);
    delete lapso;
  }else if(alg_type == 8) {
      // partition that splits all sets
      Partition part(*data);
      part.setFixed(*data,problem.get_BranchNode(),true);
      ResLP rlp(*data,part);
      Sol_Int solLB(data->getNBlock(),data->getnDestination());
      problem.get_sol(solLB);
      {
	  std::vector<double> fracSoln(data->graph.getNumNodes(),0.0);
	  std::vector<std::vector<double> > pi(data->getNPeriod(),
					       std::vector<double>(data->getnResources(),0.0));
	  rlp.solve(fracSoln,pi);
	  std::cout <<"\tRelaxed bound " << rlp.getObjVal() << " found in " << rlp.getTime()<<"\n";
      }
      rlp.solveExact(problem,solLB);
      std::cout << "\tResLP/QOL found" << solLB.obj << " in " << rlp.getTime() << " CPU sec\n";
  }
  else if(alg_type == 9){ // Window search
    // do the VLNS
    int time_limit = iterations; // an approximation for now
    // X and Y will contain the solution information when done
    obj = alg_obj->WindowSearch( disp, time_limit, show, sol, X, Y, &problem, time, shift);
    cout << "\tCompleted Window search, time points: " << time << " to " 
	 << time+shift-1 << " with objective: " << obj << endl;
    // Convert X and Y to a sol_int
    mip_sol = convertToSol_Int(X, Y, obj);
    //Sol_Int * sol_master = new Sol_Int();
    //problem.get_sol(sol_master);
    
  }   
    
  // Update the MineProblem solution
  if(obj > 0) problem.update_sol(t_max, obj, X, Y);
  //problem.current->sol_int.nT = t_max;
  //problem.current->sol_int.obj = obj;
  //problem.current->sol_int.x = X;
  //problem.current->sol_int.y = Y;

  delete sol;
  delete alg_obj;
  delete funcs;
  return 0;
}

// Update BNI with the latest solution information
// BNI is a structure within ProblemData, which sits in MineProblem
void SolverInterface::updateBNI(MineProblem &problem){

  //problem.current->sol_int.nT = t_max;
  //problem.current->sol_int.obj = t_max;
  
}

// Update a solution, required by the main parallel solver
void SolverInterface::updateSolution(Sol_Int &sol){

}

// Maybe delete solv
SolverInterface::~SolverInterface(){
  if(deleteData) delete data;
  if(mip_sol) delete mip_sol;
}


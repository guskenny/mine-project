/***************************************************************************
                            Solver.cpp
                         ----------------------------------
    last modified   : 13/5/2008
    copyright       : (C) 2013 by Dhananjay Thiruvady
    libraries	    : .
    description	    : Implementation of functions in Solver
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
#include <omp.h>

using namespace std;

#include "solver.h"

/*
  Construct the network flow class
*/

Solver::Solver(Daten &data, SolverFunctions &funcs){
  this->data = &data;
  this->funcs = &funcs;
}

/*
  Run a MIP algorithm
  Allows two different models
     In the first model a block completes at a time
     In the second model, once a jobs completes, it stays completed
*/

void Solver::MIP_run(bool simple){
  MIP *alg = new MIP(*data);
  if(simple){
    cout << "Running MIP simple ..." << endl;
    alg->MIP_simple();
  }
  else{
    cout << "Running MIP extended ..." << endl;
    alg->MIP_extended();
  }
  delete alg;
}

vector<vector<int> > Solver::createPreds(){
  vector<vector<int> > preds;
  const int nB = data->getNBlock();
  std::vector<Block> * blocks=data->getBlock();
    
  preds.resize(nB);
  for(int a=0; a<nB; a++){
    std::vector<int> * pred = (*blocks)[a].getPreds();
    int n = (*blocks)[a].getNumPred();
    preds[a].resize(n);
    for(int p=0; p<n; p++){
      int b = (*pred)[p];
      preds[a][p] = b;
    }
  }  
  
  return preds;
}

/*
  Run the network flow algorithms
  Allows three different models
*/

void Solver::NF_run(int type){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int r_max = data->getnResources();

  int penalties = false;
  // set up some basic multipliers
  vector<vector<double> > lambda;
  lambda.resize(r_max);
  for(int r = 0; r < r_max; r++){
    lambda[r].resize(t_max,0.1);
  }
  // create a generic graph
  vector<vector<int> > preds = createPreds();
  Graph *g = new Graph(nB);
  vector<vector<double> > profit;
  funcs->setLambdaProfit( profit,lambda);
  if(type == 0){
    cout << "Running network flow primal ..." << endl;
    MaxClosure_NetworkFlow_LR* net = new MaxClosure_NetworkFlow_LR(profit,*g, preds, NULL);
    vector<double> profit;
    vector<vector<double> > profit_time;
    funcs->setProfit(profit,profit_time);
    Primal(net,profit);
    delete net;
  }
  else if (type == 1){ // not implemented yet
    cout << "Running network flow dual ..." << endl;
    //alg->Dual();
  }
  else if(type == 2){ 
    cout << "Running network flow network ..." << endl;
    MaxClosure_NetworkFlow_LR* net = new MaxClosure_NetworkFlow_LR(profit,*g,preds, NULL);
    buildNetwork(net);
    delete net;
  }
  else{ 
    cout << "Running Lagrangian network flow heuristic ..." << endl;
    penalties = true;
    MaxClosure_NetworkFlow_LR* net = new MaxClosure_NetworkFlow_LR(profit,*g,preds, NULL);
    double ub = 0.0;
    bool show = true;
    runNetworkFlow(ub,show,net);
    delete net;
  }
  delete g;
}

/*
  Solve the Primal problem
  This is only to confirm that the formulations are good
  They do not use the network flow algorithm from cplex
*/

void Solver::Primal(MaxClosure_NetworkFlow_LR *net, const vector<double> &profit){
  // set up some parameters

  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  std::vector<Block> * blocks=data->getBlock();

  IloEnv env;
  IloTimer timer(env);

  set<IloNumVar> soln; // set of nodes to include
  vector<IloBoolVar> nodeVar; // variable = 1 if included in solution
  vector<IloNumVar>  edgeVar; // dual variable (flow across edge) for (pred,succ) pairs

  vector<IloBoolVar> x_i;    
  int profitSize = net->getProfitSize();

  try{
    IloModel model(env);
    // create the variables
    x_i.resize(profitSize);
    for(int i=0;i<profitSize;i++){
      x_i[i] = IloBoolVar(env,0.0,1.0,"");
    }
    
    // precedences
    IloRangeArray prec(env);
    for(int a=0; a<nB; a++){
      std::vector<int> * pred = (*blocks)[a].getPreds();
      int n = (*blocks)[a].getNumPred();
      for(int p=0; p<n; p++){
	int b = (*pred)[p];
	for(int t=0; t<t_max; t++){
	  prec.add(x_i[a] - x_i[b] <= 0);
	}
      }
    }
    model.add(prec);
    // set up the objective
    IloExpr tot(env);
    for(int i=0;i<profitSize;i++){
      tot += x_i[i] * profit[i];
    }
    model.add(IloMaximize(env,tot));
    tot.end();
    
    IloCplex cplex(model);
    cplex.setOut(env.getNullStream());
    if(!cplex.solve()){
      cout << "Executing optimization failed!" << endl;
    }

    //if(show)
    cout << "Objective value: " << cplex.getObjValue() << endl;
    
    // clean up the expressions
    prec.end();
    model.end();
  }
  catch(IloException &e){
    cout << endl << "Error during optimisation ..." << e<< endl;
  }
  catch(...){
    cout << endl << "Unknown exception during optimisation ... " << endl;
  }
}


/*
  Solve the Dual problem
  This is only to confirm that the formulations are good
  They do not use the network flow algorithm from cplex
  Need to complete this.
*/

void Solver::Dual(const vector<vector<double> > &lambda, double &upper_bound, bool show){
  // set up some parameters
}

/*
  Solve as network flow problem wwith Cplex's network algorithm
  This is the "ideal" algorithm to use
*/

void Solver::buildNetwork(MaxClosure_NetworkFlow_LR *net){
  try{
    int status;
    CPXENVptr env;
    env = CPXopenCPLEX (&status);
    if(status)
      cout << "Error opening environment ... " << status << endl;

    int status_p;
    // Create a cplex network pointer
    CPXNETptr prob = CPXNETcreateprob(env, &status_p, "NetFlowProb");
    if(status_p) cout << "Error ... " << status_p << endl;
    // set up a few things for the solver
    // make a call directly to the solver
    // create correct variables first
    int objsen = CPX_MIN;
    CPXData_LR *cd =  net->getExtNodeExtArcInfo();

    int *fromnode = &cd->from[0];
    int *tonode = &cd->to[0];
  
    CPXsetintparam(env,CPX_PARAM_SCRIND,1);
    status = CPXNETcopynet(env, prob, objsen, cd->nnodes, &cd->supply[0], NULL, cd->nedges, 
			   fromnode, tonode, &cd->low[0], &cd->up[0], &cd->obj[0], NULL);
    if(status)
      cout << "Error creating CPXNET ... " << status << endl;
    else
      cout << "Created CPXNET, now optimising ... " << endl;
    status = CPXNETprimopt(env,prob);
    if(status)
      cout << "Error during optimisation ... " << status << endl;
    else{
      double obj=0;
      CPXNETgetobjval(env,prob,&obj);
      cout << "Objective: " << obj << endl;
    }
    //cout << "Finished optimizing ... " << endl;
    //int primFeas,dualFeas;
    //CPXNETsolninfo(env,p,&primFeas,&dualFeas);
    //cout << "Primal/dual feasibility = " << primFeas << "/" << dualFeas << endl;
    //CPXNETwriteprob(env,p,"temp.net",NULL);
  }
  catch(IloException &e){
    cout << endl << "Error during optimisation ..." << e<< endl;
  }
  catch(...){
    cout << endl << "Unknown exception during optimisation ... " << endl;
  }
}

/*
  Solve as network flow problem wwith CPLEX's network algorithm
  This is the "ideal" algorithm to use
*/

vector<vector<vector<double> > > Solver::runNetworkFlow(double &upper_bound, bool show, MaxClosure_NetworkFlow_LR *net){
  // Stores the solution output
  vector<vector<vector<double> > > Yvar_sol;

  try{
    int status;
    CPXENVptr env;
    env = CPXopenCPLEX (&status);
    if(status)
      cout << "Error opening environment ... " << status << endl;

    int status_p;
    // Create a cplex network pointer
    CPXNETptr prob = CPXNETcreateprob(env, &status_p, "NetFlowProb");
    if(status_p) cout << "Error ... " << status_p << endl;
    // set up a few things for the solver
    // make a call directly to the solver
    // create correct variables first
    int objsen = CPX_MIN;
    CPXData_LR *cd =  net->getExtNodeExtArcInfo();

    int *fromnode = &cd->from[0];
    int *tonode = &cd->to[0];
  
    if(show) CPXsetintparam(env,CPX_PARAM_SCRIND,1);
    status = CPXNETcopynet(env, prob, objsen, cd->nnodes, &cd->supply[0], NULL, cd->nedges, fromnode, tonode, &cd->low[0], &cd->up[0], &cd->obj[0], NULL);
    if(status)
      cout << "Error creating CPXNET ... " << status << endl;
    else
      cout << "Created CPXNET, now optimising ... " << endl;
    status = CPXNETprimopt(env,prob);
    if(status)
      cout << "Error during optimisation ... " << status << endl;
    else{
      double obj=0;
      CPXNETgetobjval(env,prob,&obj);
      cout << "Objective: " << obj << endl;
      upper_bound = obj;
      // Get values
      vector<double> x(cd->nedges);
      status = CPXNETgetx(env, prob, &x[0], 0, cd->nedges-1);
      if(status)
	cout << "Error obtaining arc values ... " << status << endl;

      cout << "Arcs:" << endl;
      for(int i = 0; i < cd->nedges; i++){
	if(x[i] > 0.0) cout << i << ": " << x[i] << " " << endl;
      }
      vector<double> pi(cd->nnodes);
      status = CPXNETgetpi(env, prob, &pi[0], 0, cd->nnodes-1);
      if(status)
	cout << "Error obtaining node values ... " << status << endl;
      // get the values in the solution and set them to the original nodes
      cout << "Nodes:" << endl;
      int count = 0;
      for(int i = 0; i < cd->nnodes; i++){
	if(pi[i] != 0.0) {
	  count++;
	  net->setExtNodeVal(i,pi[i]);
	}
      }
      cout << "Total nodes: " << pi.size() << ", nodes used: " << count << endl;
      // Generate a solution
      Yvar_sol = funcs->generateSolution(net);
      delete cd;
    }
  }
  catch(IloException &e){
    cout << endl << "Error during optimisation ..." << e<< endl;
  }
  catch(...){
    cout << endl << "Unknown exception during optimisation ... " << endl;
  }
  return Yvar_sol;
}

vector<vector<double> > Solver::getProfit(const vector<vector<double> > &lambda){
  // initialise some parameters from the data
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const int r_max = data->getnResources();
  const double rate = data->getDiscountRate();
  std::vector<Block> * blocks=data->getBlock();

  vector<vector<double> > profit(nB);
  // Set up profit with a time component
  for(int b=0; b<nB; b++){
    profit[b].resize(t_max);
    for(int t=0; t<t_max; t++){
      profit[b][t] = 0.0;
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
      profit[b][t] = max_profit;
    }
  }
  // adjust profits
  for(int b = 0 ; b < nB ; b++) {
    for(int t=0; t<t_max-1; t++){
      profit[b][t] -= profit[b][t+1];
    }
  }

  return profit;
}

// Compute an average over a vector
double Solver::Average(const vector<double> &run_times){
  double avg = 0.0;
  for(int i=0;i<run_times.size();i++){
    avg += run_times[i];
  }
  return avg/run_times.size();
}

/*
  A single ACO iteration
  Make sure precedences are considered when building these solutions
*/

ACO_Solution* Solver::ACO(Pheromones *p, int nants, Random *r, double q_0, double lrate){
  const int nB = data->getNBlock();
  ACO_Solution *best=NULL;
  vector<ACO_Solution *> solutions(nants);
#pragma omp parallel for
  for(int k = 0; k < nants; k++){ 
    ACO_Solution *a = new ACO_Solution(this->data,funcs->profit,funcs->block_dest,funcs->succ);
    int variable = 0;   
    // vector to determine if blocks are completed or not
    vector<bool> bc(nB,false);
    // a waiting list of blocks
    vector<int> waitingBlocks;
    while(variable < nB){
      //pair<int,int> sol = 
      a->selectBlock(variable, p, r, q_0);
      variable++;
    }
    a->setBlocksWithPrec(funcs->negBlocks);
    solutions[k] = a->copy();
    delete a;
  }
  double avg = 0.0;
  for(int i = 0; i < solutions.size(); i++){
    ACO_Solution *a = solutions[i];
    //cout << a->getObj() << endl;
    avg += a->getObj();
    // Do the local pheromone update here
    /*for (int variable = 0; variable < nB; variable++){
      int time = int(a->getBlockTime(variable));
      int dest = a->getBlockDest(variable);
      if (time == -1) continue;
      p->localPheromoneUpdate(variable, time, dest, lrate);    
      }*/
    // Update the best solution
    if (i == 0){
      best = a->copy();
    }
    else if (a->getObj() > best->getObj()){
      delete best;
      best=a->copy();          
    }
    delete a;
  }
  //cout << "\t\tIB: " << best->getObj() << ", Avg: " << avg/solutions.size() << endl;
  //cout << "\t\tBest inner objective: " << best->getObj() << endl;
  // do some Local search or Beta sampling on the best solution
  /*for(int i = 0; i < 30; i++){
    ACO_Solution * a = best->betaSampling(r);
    if (a->getObj() > best->getObj()){
      delete best;
      best=a->copy();          
    }
    }*/

  return best;
}

/*
  Run an ACO
*/

ACO_Solution* Solver::ACO_run(Pheromones *p, ACO_Solution * lr_sol, int max_time, 
			      int nants, double lrate, double q_0, bool disp, 
			      bool t_limit){
  cout << "Starting ACO ..." << endl;
  int max_sols = 9;
  int gbiter = 0;
  // a random number generator
  time_t t;
  Random *r;
  r = new Random((unsigned) time(&t));
  r->next();

  ACO_Solution *gb = NULL;
  
  if (lr_sol) {
    gb = lr_sol->copy();
    if (disp) 
      cout << "\tStarting objective: " << gb->getObj() << endl;
  } 
  else{
    cout << "\tStarting without a solution ..." << endl;
    cout << "\tUsing basic heuristic" << endl;
    gb = new ACO_Solution(data,funcs->profit,funcs->block_dest,funcs->succ);
    gb->setBlocksWithPrec(funcs->negBlocks);
    cout << "\tSolution profit: " << gb->getObj() << endl;
  }
  int iteration=0;
  Timer timer;
  bool terminate = false;

  list<double> prev_ib;
  set<int> p_list;

  //cout << "\tUsing " << omp_get_num_threads() << " threads .." << endl;
  // Main loop
  while(!terminate){
    ACO_Solution *ib = ACO(p, nants, r, q_0, lrate);
    //p->display();
    //if(disp && gb) cout << "\tACO-Iter.: " << iteration << ", LB: " 
    //		<< gb->getObj() << ", CPUT: " << timer.elapsed_time(Timer::VIRTUAL) 
    //			<< ", RT: " << timer.elapsed_time(Timer::REAL) << endl;
    if(gb == NULL){ // in the first iteration we initialize all the variables
      gb = ib->copy();
      gbiter=iteration;
    }
    else{ 
      if(prev_ib.size()>max_sols) prev_ib.pop_front();
      prev_ib.push_back(ib->getObj());
      if(ib->getObj() > gb->getObj()){
	delete gb;
	gb=ib->copy();
	gbiter=iteration;
	cout << "\tACO-Iter.: " << iteration << ", LB: " 
	     << gb->getObj() << ", CPUT: " << timer.elapsed_time(Timer::VIRTUAL) 
	     << ", RT: " << timer.elapsed_time(Timer::REAL) << endl;
	
      }
      p->globalPheromoneUpdate(gb->getBlockTimes(), gb->getBlockDests(), lrate, gb->getObj());
      double tempo = prev_ib.front();
      bool same = true;
      if(prev_ib.size() > max_sols){ // check for convergence, if 5 solutions are the same, assume converged
	for (list<double>::iterator sit = prev_ib.begin(); sit != prev_ib.end(); sit++) {
	  if((*sit)!=tempo) {
	    same =false;
	    break;
	  }
	}
	if(same) { // Reset the pheromones
	  //cout << "Pheromones converged ... re-setting" << endl; 
	  p->reset(this->data);
	  prev_ib.clear();
	}
      }
    }
    iteration++;
    if (!t_limit && iteration > max_time) terminate = true;
    else if( t_limit && timer.elapsed_time(Timer::VIRTUAL) > max_time) terminate = true;
    delete ib;
  }
  if(gb){
    cout << "\tBest objective: " << gb->getObj() << endl;
  }
  delete r;
  return gb;
}

/*
  Run an ACO will an interface to a MineProblem
*/

ACO_Solution* Solver::ACO_run(Pheromones *p, ACO_Solution * lr_sol, int max_time, 
			      int nants, double lrate, double q_0, bool disp, 
			      bool t_limit, MineProblem *problem){ 
  cout << "\tStarting ACO ..." << endl;
  int max_sols = 9;
  int gbiter=0;
  // a random number generator
  time_t t;
  Random *r;
  r = new Random((unsigned) time(&t));
  r->next();

  ACO_Solution *gb = NULL;
  
  if (lr_sol) {
    gb = lr_sol->copy();
    if (disp) 
      cout << "\tStarting objective: " << gb->getObj() << endl;
  } 
  else{
    cout << "\tStarting without a solution ..." << endl;
    cout << "\tUsing basic heuristic" << endl;
    gb = new ACO_Solution(data,funcs->profit,funcs->block_dest,funcs->succ);
    gb->setBlocksWithPrec(funcs->negBlocks);
    cout << "\tSolution profit: " << gb->getObj() << endl;
  }
  int iteration=0;
  Timer timer;
  bool terminate = false;

  list<double> prev_ib;
  set<int> p_list;

  //cout << "\tUsing " << omp_get_num_threads() << " threads .." << endl;
  // Main loop
  while(!terminate){
    ACO_Solution *ib = ACO(p, nants, r, q_0, lrate);
    //p->display();
    //if(disp && gb) cout << "\tACO-Iter.: " << iteration << ", LB: " 
    //			<< gb->getObj() << ", CPUT: " << timer.elapsed_time(Timer::VIRTUAL) 
    //			<< ", RT: " << timer.elapsed_time(Timer::REAL) << endl;
    if(gb == NULL){ // in the first iteration we initialize all the variables
      gb = ib->copy();
      gbiter=iteration;
    }
    else{ 
      if(prev_ib.size()>max_sols) prev_ib.pop_front();
      prev_ib.push_back(ib->getObj());
      if(ib->getObj() > gb->getObj()){
	delete gb;
	gb=ib->copy();
	gbiter=iteration;
	cout << "\tACO-Iter.: " << iteration << ", LB: " 
	     << gb->getObj() << ", CPUT: " << timer.elapsed_time(Timer::VIRTUAL) 
	     << ", RT: " << timer.elapsed_time(Timer::REAL) << endl;
	// update the problem solution if we have found something better
	if(problem != NULL){
	  Sol_Int *sol = new Sol_Int();
	  funcs->generateParSolution(ib, *sol);
	  problem->update_sol(*sol);
	  problem->set_updated();
	  delete sol;
	}
	/*vector<int> X;// = problem.current->sol_int.x;
	vector<vector<double> > Y;// = problem.current->sol_int.y;
	double obj;
	problem.get_sol(time, obj, X, Y);
	cout << "\t\tNew updated ACO .." << obj << endl;*/
      }
      p->globalPheromoneUpdate(gb->getBlockTimes(), gb->getBlockDests(), lrate, gb->getObj());
      double tempo = prev_ib.front();
      bool same = true;
      if(prev_ib.size() > max_sols){ // check for convergence, if 5 solutions are the same, assume converged
	for (list<double>::iterator sit = prev_ib.begin(); sit != prev_ib.end(); sit++) {
	  if((*sit)!=tempo) {
	    same =false;
	    break;
	  }
	}
	if(same) { // Reset the pheromones
	  cout << "\tPheromones converged ... re-setting" << endl; 
	  p->reset(this->data);
	  prev_ib.clear();
	}
      }
    }
    // check to see if we have a new global bes
    if(problem != NULL){
      if(problem->get_updated()){ // yes we do
	vector<int> X;// = problem.current->sol_int.x;
	vector<vector<double> > Y;// = problem.current->sol_int.y;
	double obj;
	int time;
	problem->get_sol(time, obj, X, Y);
	// Construct a solution from problem or X and Y 
	ACO_Solution* sol = new ACO_Solution(data,funcs->profit,
					     funcs->block_dest,funcs->succ,
					     funcs->negBlocks,
					     X, Y, false);
	//cout << "Solution objective : " << sol->getObj() << endl;
	if(gb->getObj() < sol->getObj()){
	  cout << "\t\tUpdating best ACO solution ... previous best: " 
	       << gb->getObj() << ", new obj: " << sol->getObj() << endl;
	  delete gb;
	  gb = sol->copy();
	}
	delete sol;
	//cout << "\tUpdated solution on the ACO side." << endl;
      }
    }

    iteration++;
    if (!t_limit && iteration > max_time) terminate = true;
    else if( t_limit && timer.elapsed_time(Timer::VIRTUAL) > max_time) terminate = true;
    delete ib;
  }
  if(gb){
    cout << "\tBest objective: " << gb->getObj() << endl;
  }
  delete r;
  return gb;
}


/*
  A Lagrangian relaxation algorithm
  A basic implementation at this stage, no parametric optimisation has been investigated
  
*/

ACO_Solution* Solver::Lagrangian_heuristic(bool disp, double time_limit, bool show, 
					   const double &bound, int iterations, int nants, double lrate, 
					   double q_0, bool use_aco, MineProblem *problem){
  // initialise some parameters from the data
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int r_max = data->getnResources();

  // Could use a MIP or another type of algorithm for the relaxation
  //MIP *alg = new MIP(this->data);
  //MaxClosure_NetworkFlow_LRFlow *nf_alg = new MaxClosure_NetworkFlow_LRFlow(this->data);

  int itr = 0;
  double  best_upper_bound = 1e99;
  double  best_lower_bound = -1e99;
  double gap = 1.0;
  double ub = 1e99;
  double delta = 2.0;
  int ubc =0;
  const int max_iterations = 1000; 
  //const int max_iterations = 2; 

  vector<vector<double> > lambda;
  vector<vector<double> > best_lambda; 

  if(disp)  cout << "Starting Lagrangian relaxation heuristic ... " << endl;

  funcs->setLambda(lambda);
  best_lambda = lambda;
  // Create a network object
  // Usable by any algorithm
  vector<vector<int> > preds = createPreds();
  Graph *g = new Graph(nB);
  vector<vector<double> > profit;
  funcs->setLambdaProfit(profit,lambda);
  MaxClosure_NetworkFlow_LR* net = new MaxClosure_NetworkFlow_LR(profit,*g,preds,NULL);
  // To implement timer related stuff
  Timer timer;
  vector<double> run_times;

  bool terminate = false;
  // some parameters for updating delta
  int last_ub_change = 0;
  bool increase = true;

  // To record the best solution
  ACO_Solution *best_sol = new ACO_Solution(data,funcs->profit,funcs->block_dest,funcs->succ);

  // set up a pheromone matrix to be learned through the execution of the algorithm
  Pheromones *p = new Pheromones(this->data, funcs->profit);

  // main loop
  while(itr < max_iterations && delta > 0.01 && gap > 0.01 && !terminate){
    ub=0.0;
    //solve the Langrangian dual problem for every machine, line 2
    double upper_bound=0.0;
    //vector<vector<vector<double> > > Yvar = alg->MIP_relaxed(lambda, upper_bound, show);
    //vector<vector<vector<double> > > Yvar = alg->MIP_relaxed_efficient(lambda, upper_bound, show);
    //vector<vector<vector<double> > > Yvar = nf_alg->LagrangianRelaxationNetwork(upper_bound, show, net);
    vector<vector<vector<double> > > Yvar;// = nf_alg->Test(upper_bound, show, net);
    int status = 0;
    double real_time_taken = timer.elapsed_time(Timer::REAL);
    if(itr > 0) {// Update the profits in the network with the new lambdas
      vector<vector<double> > profit = getProfit(lambda);
      net->updateProfit(profit);
      // now set up the supply's from the updated information
      status = net->updateSupply();
      if(status)
	cout << "Error when updating supplies ... " << status << endl;
    }
    status = net->solve();
    double actual_run_time = timer.elapsed_time(Timer::REAL) - real_time_taken; 
    run_times.push_back(actual_run_time);

    if(status)
      cout << "Error during optimisation ... " << status << endl;
    else{      // Generate a solution
      upper_bound = net->getSolInfo();
      Yvar = funcs->generateSolution(net);
    }
    // Make sure the net is clear of any data from the MIP
    net->resetValues();

    // Determine the objective after adding the constant
    double sum = 0.0;
    for(int r = 0; r < r_max ; r++){
      for(int t = 0; t < t_max ; t++){
	sum += data->getLimit(r, t) * lambda[r][t];
      }
    }
    upper_bound = (upper_bound*1e8+sum);
    cout << "\tLR objective (with constant): " <<  upper_bound << endl;

    // Run an ACO instead
    ACO_Solution * lr_sol = new ACO_Solution(data,funcs->profit,funcs->block_dest,funcs->succ);
    lr_sol->setTimes(Yvar);
    int shift = 0;
    double obj = lr_sol->setBlocks(funcs->negBlocks, shift);
    //double obj = lr_sol->setBlocksPar();
    cout << "\tRepaired objective: " << obj << endl;
    double lb = lr_sol->getObj(); 
    ACO_Solution *ib = NULL;
    // update the best solution
    if(best_sol){
      if(lr_sol->getObj() > best_sol->getObj()){
	delete best_sol;
	best_sol = lr_sol->copy();
      }
    }
    else {
      best_sol = lr_sol->copy();
    }
    if(use_aco){
      ib = ACO_run(p, lr_sol, iterations, nants, lrate, q_0, disp, false, problem);
      lb = ib->getObj(); 
    }
    else{
      ib = lr_sol->copy();
    }
    
    // update the best solution
    if(best_sol){
      if(ib->getObj() > best_sol->getObj()){
	delete best_sol;
	best_sol = ib->copy();
	// update the problem solution if we have found something better
	Sol_Int *sol = new Sol_Int();
	funcs->generateParSolution(lr_sol, *sol);
	problem->update_sol(*sol);
	problem->set_updated();
	//cout << "\t\tNew ACO best found: updating master .." << endl;
	delete sol;
	cout << "\tBest solution updated, iteration " << itr << endl;
      }
    }
    else {
      cout << "\tBest solution updated, iteration " << itr << endl;
      best_sol = ib->copy();
    }

    // check to see if we have a new best solution from any other thread
    if(problem != NULL){
      if(problem->get_updated()){ // yes we do
	vector<int> X;// = problem.current->sol_int.x;
	vector<vector<double> > Y;// = problem.current->sol_int.y;
	double obj;
	int time;
	problem->get_sol(time, obj, X, Y);
	// Construct a solution from problem or X and Y 
	ACO_Solution* sol = new ACO_Solution(data,funcs->profit,
					     funcs->block_dest,funcs->succ,
					     funcs->negBlocks,
					     X, Y, false);
	if(best_sol->getObj() < sol->getObj()){
	  cout << "\t\tUpdating best LR solution..." << endl;
	  delete best_sol;
	  best_sol = sol->copy();
	}
	delete sol;
	problem->set_updated();
      }
    }

    //if(bound != -1) lb = bound;
    //cout << "Lower bound " << lb << endl;
    if(lb > best_lower_bound) best_lower_bound = lb;
    //double old_upper_bound = upper_bound;
    //if(disp) cout << "\n\t" << "Cost of relaxed problem" << ": " << upper_bound << endl; 
    ub += upper_bound;

    // Decrease delta every fifth iteration
    // Increase delta if the solution has been improving
    ubc++;
    if(ubc>=5){
      ubc=0;
      delta*=0.95;
      increase = false;
    }
    else if(itr - last_ub_change < 5 && !increase){
      increase = true;
      delta += delta + 0.1*delta;
    }

    if( itr == 0 ) {
      best_lambda = lambda;
      best_upper_bound = ub;
    }
    else if(best_upper_bound > ub){
      ubc=0;
      best_lambda = lambda;
      best_upper_bound = ub;
      last_ub_change = itr;
    }
    // update the multipliers
    double relaxed_gap =  (best_upper_bound-best_lower_bound);
    lambda = funcs->update_lambda(relaxed_gap,lambda,delta,Yvar);
    if(status) 
      cout << "Problem updating supplies ... " << endl;
    // recompute the gap
    gap = (double)(best_upper_bound-best_lower_bound)/best_upper_bound;
    
    if(disp) {
      cout << setprecision(12) << "\n\tIter. " << itr+1 << ": " 
	   << "Gap: " << gap << ", LB: " <<  best_lower_bound 
	   << ", UB: " << best_upper_bound << ", RUB: " << ub 
	   << ", Delta: " << delta;  

      cout << setprecision(5) 
	   << " , CPUT: " << timer.elapsed_time(Timer::VIRTUAL)
	   << " , RT: " << timer.elapsed_time(Timer::REAL) 
	   << " , ART: " << actual_run_time 
	   <<  endl; 
    }
    delete ib;
    delete lr_sol;
    itr++;
    if(timer.elapsed_time(Timer::REAL) > time_limit) terminate = true;
  } 
  
  cout << gap << "\t" << best_lower_bound << "\t" << best_upper_bound 
       << "\t" << timer.elapsed_time(Timer::VIRTUAL) 
       << "\t" << timer.elapsed_time(Timer::REAL) 
       << "\t" << Average(run_times) << endl;
  
  delete g;
  delete net;
  delete p;
  //delete alg;
  //delete nf_alg;
  lambda.clear();
  best_lambda.clear();
  return best_sol;
}

/*
  A large neighbourhood search
  Given a set of fixed variables, run a MIP.
*/

double Solver::LargeNeighbourhoodSearch( bool disp, 
					 int time_limit, bool show, ACO_Solution *sol,
					 vector<int> &X, 
					 vector<vector<double> > &Y,
					 MineProblem *problem
					 ){
  cout << "\tCommencing a very large neighbourhood search ... " << endl; 
  // run an ACO and get its solution
  double iterations = 10;
  int nants = 10;
  double lrate = 0.2;
  double q_0 = 0.1;
  bool use_tl = false;
  Pheromones *p = new Pheromones(this->data, funcs->profit);
  ACO_Solution * a = ACO_run(p, sol, iterations, nants, lrate, q_0, disp, use_tl, problem);
  vector<vector<vector<double> > > Y_var = a->mappedYvar();
  delete a;
  delete p;

  Timer timer;
  double time_taken = timer.elapsed_time(Timer::VIRTUAL);
  double r_time_taken = timer.elapsed_time(Timer::REAL);
  MIP *alg = new MIP(*data);
  // Due to model considerations, call a function in mip.cpp
  //double best_obj = alg->VLNS(Y_var, funcs->profit, timer, problem);
  //double best_obj = alg->VLNS_BW(Y_var, funcs->profit, timer, problem);
  //cout << "Commencing VLNS ... " << endl;
  double best_obj = alg->VLNS_Small(Y_var, funcs->profit, timer, X, Y, problem);
  // report the best solution found
  cout << "\tCompleted VLNS: " << best_obj << ", CPUT: " 
       << timer.elapsed_time(Timer::VIRTUAL) - time_taken << ", RT: "  
       << timer.elapsed_time(Timer::REAL) - r_time_taken
       << endl;
  delete alg;
  return best_obj;
}


/*
  A large neighbourhood search
  Given a set of fixed variables, run a MIP.
*/

double Solver::WindowSearch( bool disp, 
			     int time_limit, bool show, ACO_Solution *sol,
			     vector<int> &X, 
			     vector<vector<double> > &Y,
			     MineProblem *problem,
			     int time, int shift
			     ){
  cout << "\tCommencing a window search ... time: " << time << ", shift: " << shift << endl; 
  // run an ACO and get its solution
  double iterations = 10;
  int nants = 10;
  double lrate = 0.2;
  double q_0 = 0.1;
  bool use_tl = false;
  Pheromones *p = new Pheromones(this->data, funcs->profit);
  ACO_Solution * a = ACO_run(p, sol, iterations, nants, lrate, q_0, disp, use_tl, problem);
  vector<vector<vector<double> > > Y_var = a->mappedYvar();
  delete a;
  delete p;

  Timer timer;
  double time_taken = timer.elapsed_time(Timer::VIRTUAL);
  double r_time_taken = timer.elapsed_time(Timer::REAL);
  MIP *alg = new MIP(*data);
  double best_obj = alg->VLNS_Window(Y_var, funcs->profit, timer, X, Y, problem, time, shift);
  // report the best solution found
  cout << "\tCompleted VLNS window: " << best_obj << ", CPUT: " 
       << timer.elapsed_time(Timer::VIRTUAL) - time_taken << ", RT: "  
       << timer.elapsed_time(Timer::REAL) - r_time_taken
       << endl;
  delete alg;
  return best_obj;
}

/*
  Fix some variables when doing the large neighbouhood search
 */

double Solver::LNS_fixed(bool disp, 
			 int time_limit, bool show,
			 MineProblem *problem,
			 vector<int> &X, 
			 vector<vector<double> > &Y 
			 ){
  cout << "\tCommencing a very large neighbourhood search (with fixing variables) ... " << endl; 
  Timer timer;
  double time_taken = timer.elapsed_time(Timer::VIRTUAL);
  double r_time_taken = timer.elapsed_time(Timer::REAL);
  MIP *alg = new MIP(*data);
  // Due to model considerations, call a function in mip.cpp
  double best_obj = alg->VLNS_fixed(funcs->profit, timer, X, Y, problem);
  // report the best solution found
  cout << "\tCompleted VLNS: " << best_obj << ", CPUT: " 
       << timer.elapsed_time(Timer::VIRTUAL) - time_taken << ", RT: "  
       << timer.elapsed_time(Timer::REAL) - r_time_taken
       << endl;
  delete alg;
  return best_obj;
}

/*
  Run an algorithm specified by type.
*/

bool Solver::run_algorithm(bool disp, int time_limit, bool show, int type, 
			   bool simple, int nf_type, double bound, int nants, double lrate, 
			   double q_0, bool use_aco){
  if(type == 1){ // run a MIP
    MIP_run(simple);
  }
  else if (type == 2){ // run a Lagrangian relaxation
    // Use the time limit as iterations
    int iterations = time_limit;
    // Check if a mistake has been made in the input
    //   Is easy enough to make a mistake since times are in seconds
    if(iterations > 100) iterations = 100;
    MineProblem *problem = NULL;
    // fix some dummy node info, i.e., every task may only be done at the first time_point
    /*for(int b = 0; b < nB; b++){
      std::array<int, 2> tim = {0,1};
      problem->current->node.time.push_back(tim);
      // the statements below are not needed, just testing
      problem->current->node.time[b][0]=0;
      problem->current->node.time[b][1]=t_max;
      }*/
    //cout << "Finished determining BNI ... " << endl;
    ACO_Solution *a = Lagrangian_heuristic(disp, time_limit, show, bound, iterations, nants, lrate, q_0, use_aco, problem);
    delete a;
  }
  else if (type == 3){ // run a network flow
    NF_run(nf_type);  
  }
  else if (type == 4){ // run an ACO
    // Use the time limit as iterations
    int iterations = time_limit;
    // Check if a mistake has been made in the input
    //   Is easy enough to make a mistake since times are in seconds
    if(iterations > 5000) iterations = 5000;
    Pheromones *p = new Pheromones(this->data, funcs->profit);
    ACO_Solution * a = ACO_run(p, NULL, iterations, nants, lrate, q_0, disp, false);
    delete a;
    delete p;
  }
  else if (type == 5){ // LR and ACO in parallel
    // Use the time limit as iterations
    int iterations = time_limit;
    // Check if a mistake has been made in the input
    //   Is easy enough to make a mistake since times are in seconds
    if(iterations > 5000) iterations = 5000;
    Pheromones *p = new Pheromones(this->data, funcs->profit);
    ACO_Solution * a = ACO_run(p, NULL, iterations, nants, lrate, q_0, disp, false);
    delete a;
    delete p;
  }
  else if (type == 6){ // LR and ACO in parallel
    // X and Y will store the results
    vector<int> X;
    vector<vector<double> > Y;
    MineProblem *problem = NULL;
    LargeNeighbourhoodSearch(disp, time_limit, show, NULL, X, Y, problem);
  }
  return true;
  
}

Solver::~Solver(){
}

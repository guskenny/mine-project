/*************************************************************************** 
                          Implementation of MIP algorithms 
                         ---------------------------------- 
    last modified   : 13/5/2016 
    copyright       : (C) 2016 by Dhananjay Thiruvady 
    libraries		: . 
    description		: various MIP algorithms
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
#include <vector> 
#include <cmath> 
#include <string.h>
#include <sstream>
#include <iomanip>
 
using namespace std; 

#include "mip.h" 
#include "parameters.h" 

// set up data in the constructor
MIP::MIP(Daten &dat){
  this->data = &dat;
}

/*
  Set up the variables
  This includes X and Y
*/

void MIP::setVariables( GRBModel &model, 
			vector<vector<GRBVar> > &Xvar,  
			vector<vector<vector<GRBVar> > > &Yvar){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  // Create variables and set them to be binary
  for(int b=0; b<nB; b++){
    Xvar[b].resize(t_max);
    Yvar[b].resize(t_max);
    for(int t=0; t<t_max; t++){
      Xvar[b][t] = model.addVar(0,1,0,GRB_BINARY);
      Yvar[b][t].resize(d_max);
      for(int d=0; d<d_max; d++){
	Yvar[b][t][d] = model.addVar(0,1,0,GRB_CONTINUOUS);
      }
    }
  }
  // Update model to integrate the variables
  model.update();
}

/*
  Set the objective
*/

void MIP::setObjective( GRBModel &model, 
			vector<vector<vector<GRBVar> > > &Yvar){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const double rate = data->getDiscountRate();
  std::vector<Block> * blocks=data->getBlock();
  // objective function
  cout << endl << "Setting up the objective ..." << endl;
  GRBLinExpr exprObj=0;
  for(int b=0; b<nB; b++){
    for(int d=0; d<d_max; d++){
      double coef = (*blocks)[b].getProfit(d);
      for(int t=0; t<t_max; t++){
	exprObj += coef*Yvar[b][t][d];
	  coef /= (1+rate); 
      }
    }
  }
  //set the objective and optimise
  model.setObjective(exprObj,GRB_MAXIMIZE);
  model.update();
}

/*
  Set up the precedence constraints
*/

void MIP::precedenceConstraints( GRBModel &model, 
				 vector<vector<GRBVar> > &Xvar){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  std::vector<Block> * blocks=data->getBlock();

  ///Precedence (7) 
  for(int a=0; a<nB; a++){
    std::vector<int> * pred = (*blocks)[a].getPreds();
    int n = (*blocks)[a].getNumPred();
    for(int p=0; p<n; p++){
      int b = (*pred)[p];
      for(int t=0; t<t_max; t++){
	GRBLinExpr exprLhs = 0;
	GRBLinExpr exprRhs = 0;
	for(int j=0; j<=t; j++){
	  exprLhs += Xvar[a][j];
	  exprRhs += Xvar[b][j];
	}
	model.addConstr(exprLhs <= exprRhs, "");
      }
    }
  }
  model.update();

}

/*
  Set up the precedence constraints for the extended model
*/

void MIP::precedenceConstraintsExt( GRBModel &model, 
				    vector<vector<GRBVar> > &Xvar){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  std::vector<Block> * blocks=data->getBlock();
  
  ///Precedence (7) 
  for(int a=0; a<nB; a++){
    std::vector<int> * pred = (*blocks)[a].getPreds();
    int n = (*blocks)[a].getNumPred();
    for(int p=0; p<n; p++){
      int b = (*pred)[p];
      for(int t=0; t<t_max; t++){
	model.addConstr(Xvar[a][t] <= Xvar[b][t], "");
      }
    }
  }
}

/*
  Sum Destination and clique block constraints

*/

void MIP::sumDestCliqueBlockConstraint( GRBModel &model, 
					vector<vector<GRBVar> > &Xvar,
					vector<vector<vector<GRBVar> > > &Yvar){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();

  // SumDest (8) 
  for(int b=0; b<nB; b++){
    for(int t=0; t<t_max; t++){
      GRBLinExpr expr = 0;
      expr -= Xvar[b][t];
      for(int d=0; d<d_max; d++)
	expr += Yvar[b][t][d];
      model.addConstr(expr == 0, "");
    }
  }
  model.update();

  // CliqueBlock (9) 
  for(int b=0; b<nB; b++){
    GRBLinExpr expr = 0;
    for(int t=0; t<t_max; t++)
      expr += Xvar[b][t];
    model.addConstr(expr <= 1, "");
  }
  model.update();
}

/*
 Sum destination and block completion constraints
 For the extended model
*/

void MIP::sumDestCompletionConstraint( GRBModel &model, 
				       vector<vector<GRBVar> > &Xvar,
				       vector<vector<vector<GRBVar> > > &Yvar){
  
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();

  // SumDest (8) 
  for(int b=0; b<nB; b++){
    for(int t=1; t<t_max; t++){
      GRBLinExpr expr = 0;
      expr = Xvar[b][t]-Xvar[b][t-1];
      GRBLinExpr expr2 = 0;
      for(int d=0; d<d_max; d++)
	expr2 += Yvar[b][t][d];
      model.addConstr(expr2 == expr, "");
    }
    GRBLinExpr expr2 = 0;
    for(int d=0; d<d_max; d++)
      expr2 += Yvar[b][0][d];
    model.addConstr(expr2 == Xvar[b][0], "");
  }
  model.update();
  
  // changed model, so once a job is completed it stays that way
  for(int b=0; b<nB; b++){
    for(int t=0; t<t_max-1; t++){
      model.addConstr(Xvar[b][t] <= Xvar[b][t+1], "");
    }
  }
  model.update();
}

/*
  Resource constraints
*/

void MIP::resourceConstraints( GRBModel &model, 
			       vector<vector<vector<GRBVar> > > &Yvar){

  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const int r_max = data->getnResources();
  std::vector<Block> * blocks=data->getBlock();

  // Resource constraints (10) 
  for(int r=0; r<r_max; r++){
    for(int t=0; t<t_max; t++){
      char cType = data->getResConstrType(r, t);
      GRBLinExpr expr = 0;
      for(int b=0; b<nB; b++){
	for(int d=0; d<d_max; d++){
	  double coef = (*blocks)[b].getRCoef(d,r);
	  expr += coef*Yvar[b][t][d];
	}
      }
      if(cType == 'L'){ 
	model.addConstr(expr <= data->getLimit(r, t), "");
	model.addConstr(-GRB_INFINITY <= expr, "");
      }else if(cType == 'R'){
	model.addConstr(expr >= data->getLimit(r, t), "");
	model.addConstr(expr <= GRB_INFINITY, "");
      }else if(cType == 'I'){
	model.addConstr(expr <= data->getLimit(r, t, 1), "");
	model.addConstr(expr >= data->getLimit(r,t), "");
      }
    }
  }
  model.update();
}


/*
  Solve a MIP (simple model) for the OPBS
  Data is read from the implementation by Davaa: Daten
  
  This implementation is the simlpe model where x_bt = 1 
    when a block is mined at time point t. It is less effecient 
    than the second model where, once a block is mined, it 
    stays mined.
 */

void MIP::MIP_simple(){
  // set up some parameters
  const int nB = data->getNBlock();
  // set up the model now
  GRBEnv *env = NULL;
  vector<vector<GRBVar> > Xvar(nB);
  vector<vector<vector<GRBVar> > > Yvar(nB);
  
  cout << "Running the simple model ... " << endl;

  try {
    env = new GRBEnv();
    GRBModel model = GRBModel(*env);
    //model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_StringAttr_ModelName, "MIP_OPBS");
    setVariables(model,Xvar,Yvar);
    // Constraints 
    precedenceConstraints(model,Xvar);
    sumDestCliqueBlockConstraint(model,Xvar,Yvar);
    resourceConstraints(model,Yvar);
    // Set up the objective
    setObjective(model,Yvar);
    // optimize
    model.optimize();
    // display the solution
    cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
  } 
  catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } 
  catch(...) {
    cout << "Exception during optimization" << endl;
  }
  if(env)
    delete env;
}

/*
  Solve a MIP (extended model) for the OPBS
  Data is read from the implementation by Davaa: Daten
  
  This implementation is the extended model where x_bt = 1 
    when a block is mined at time point t and stays mined. 
    It is more effecient than the earlier model as evidenced
    in multiple studies, e.g., project scheduling, resource
    constrained job scheduling, etc.
 */

void MIP::MIP_extended(){
  // set up some parameters
  const int nB = data->getNBlock();

  // set up the model now
  GRBEnv *env = NULL;
  vector<vector<GRBVar> > Xvar(nB);
  vector<vector<vector<GRBVar> > > Yvar(nB);
  
  cout << endl << "Running the extended model ... " << endl;
  try {
    env = new GRBEnv();
    GRBModel model = GRBModel(*env);
    //model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_StringAttr_ModelName, "MIP_OPBS");
    setVariables(model,Xvar,Yvar);
    // Constraints
    precedenceConstraintsExt(model,Xvar);
    sumDestCompletionConstraint(model,Xvar,Yvar);
    resourceConstraints(model,Yvar);
    // Set up the objective
    setObjective(model,Yvar);
    model.optimize();
    // display the solution
    cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
  } 
  catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } 
  catch(...) {
    cout << "Exception during optimization" << endl;
  }
  if (env) 
    delete env;
}

/*
  Solve a MIP (relaxed model) for the OPBS
  Data is read from the implementation by Davaa: Daten
  
  This implementation is the extended model where x_bt = 1 
    when a block is mined at time point t and stays mined. 
    It is more effecient than the earlier model as evidenced
    in multiple studies, e.g., project scheduling, resource
    constrained job scheduling, etc.

  We solve the relaxed problem here with Lagrangian multipliers
    Lambda. The resource constraint is relaxed and added as a 
    term to the objective.
 */

vector<vector<vector<double> > > MIP::MIP_relaxed(const vector<vector<double> > &lambda, double &upper_bound, bool show){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const int r_max = data->getnResources();
  const double rate = data->getDiscountRate();
  std::vector<Block> * blocks=data->getBlock();

  // set up the model now
  GRBEnv *env = NULL;
  // the variables
  vector<vector<GRBVar> > Xvar(nB);
  vector<vector<vector<GRBVar> > > Yvar(nB);
  // the solution to be returned, output of Yvar
  vector<vector<vector<double> > > Yvar_sol(nB);
  
  try {
    env = new GRBEnv();
    GRBModel model = GRBModel(*env);
    if(!show) model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_StringAttr_ModelName, "MIP_OPBS");
    setVariables(model,Xvar,Yvar);
    // Set up the solution
    for(int b=0; b<nB; b++){
      Yvar_sol[b].resize(t_max);
      for(int t=0; t<t_max; t++){
	Yvar_sol[b][t].resize(d_max,0.0);
      }
    }
    //Constraints
    precedenceConstraintsExt(model, Xvar);
    sumDestCompletionConstraint( model, Xvar, Yvar); 
    // the next component is to include the penalty into the costs
    GRBLinExpr exprObjLag=0;
    for(int r=0; r<r_max; r++){
      for(int t=0; t<t_max; t++){
	char cType = data->getResConstrType(r, t);
	GRBLinExpr expr = 0;
	for(int b=0; b<nB; b++){
	  for(int d=0; d<d_max; d++){
	    double coef = (*blocks)[b].getRCoef(d,r);
	    expr -= coef*Yvar[b][t][d];
	  }
	}
	if(cType == 'L'){ 
	  expr += data->getLimit(r, t);
	  exprObjLag += lambda[r][t] * expr; 
	}
	else if(cType == 'R'){ // not really used for our datasetd
	  exprObjLag += lambda[r][t] * (data->getLimit(r, t) - expr); 
	}
	else if(cType == 'I'){// not really used for our datasetd
	  expr -= data->getLimit(r, t, 1);
	  exprObjLag += lambda[r][t] * (data->getLimit(r, t) - expr); 
	}
      }
    }

    // objective function
    cout << endl << "Setting up the objective ..." << endl;
    GRBLinExpr exprObj=0;
    for(int b=0; b<nB; b++){
      for(int d=0; d<d_max; d++){
	double coef = (*blocks)[b].getProfit(d);
	for(int t=0; t<t_max; t++){
	  exprObj += coef*Yvar[b][t][d];
	  coef /= (1+rate); 
	}
      }
    }
    // include the penalties
    exprObj += exprObjLag;
    //set the objective
    model.setObjective(exprObj,GRB_MAXIMIZE);
    model.update();
    // optimize
    model.optimize();
    // display the solution
    if(show) cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
    upper_bound = model.get(GRB_DoubleAttr_ObjVal);

    // output solution and also capture the Y variables    
    for(int b=0; b<nB; b++){
      bool YesNo = 0;
      for(int t =0; t<t_max; t++){
	if(Xvar[b][t].get(GRB_DoubleAttr_X)>0.9){
	  for(int d =0; d<d_max; d++){
	    Yvar_sol[b][t][d] = Yvar[b][t][d].get(GRB_DoubleAttr_X);
	  }
	  YesNo = 1;
	  break;
	}
      }
    }
  } 
  catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } 
  catch(...) {
    cout << "Exception during optimization" << endl;
  }
  if(env)
    delete env;
  return Yvar_sol;
}


/*
  Solve a MIP (relaxed and efficient model) for the OPBS
  Data is read from the implementation by Davaa: Daten
  
  This implementation is the extended model where x_bt = 1 
    when a block is mined at time point t and stays mined. 
    It is more effecient than the earlier model as evidenced
    in multiple studies, e.g., project scheduling, resource
    constrained job scheduling, etc.

  We solve the relaxed problem (only with x variables) here 
    with Lagrangian multipliers Lambda. The resource constraint 
    is relaxed and added as a term to the objective.
 */

vector<vector<vector<double> > > MIP::MIP_relaxed_efficient(const vector<vector<double> > &lambda, double &upper_bound, bool show){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const int r_max = data->getnResources();
  const double rate = data->getDiscountRate();
  std::vector<Block> * blocks=data->getBlock();

  // set up the model now
  GRBEnv *env = NULL;
  // the variables
  vector<vector<GRBVar> > Xvar(nB);
  // the solution to be returned, output of Yvar
  vector<vector<vector<double> > > Yvar_sol(nB);

  // set up a profit vector for th objective
  vector<vector<double> > profit(nB);
  vector<vector<int> > max_block_profit(nB); // destination if block b in time t
  for(int b=0; b<nB; b++){
    profit[b].resize(t_max,0.0);
    max_block_profit[b].resize(t_max,0);
    for(int t=0; t<t_max; t++){
      // first determine the desination with the max profit
      double max_profit = -1e99;
      int dest = 0;
      for(int d=0; d<d_max; d++){
	double t_profit = (*blocks)[b].getProfit(d)/pow(1+rate,t); 
	for(int r=0; r<r_max; r++)
	  t_profit -= (lambda[r][t] * (*blocks)[b].getRCoef(d,r));
	if(t_profit > max_profit) {
	   max_profit = t_profit;
	   dest = d;
	}
      }
      profit[b][t] = max_profit;
      max_block_profit[b][t] = dest;
    }
  }
  
  try {
    env = new GRBEnv();
    GRBModel model = GRBModel(*env);
    if(!show) model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_StringAttr_ModelName, "MIP_OPBS");
    // Create variables and set them to be binary
    for(int b=0; b<nB; b++){
      Xvar[b].resize(t_max);
      Yvar_sol[b].resize(t_max);
      for(int t=0; t<t_max; t++){
	Xvar[b][t] = model.addVar(0,1,0,GRB_BINARY);
	Yvar_sol[b][t].resize(d_max,0.0);
      }
    }
    // Update model to integrate the variables
    model.update();

    // Constraints
    precedenceConstraintsExt(model,Xvar);
    
    // changed model, so once a job is completed it stays that way
    for(int b=0; b<nB; b++){
      for(int t=0; t<t_max-1; t++){
	model.addConstr(Xvar[b][t] <= Xvar[b][t+1], "");
      }
    }
    model.update();

    // objective function
    if(show) cout << endl << "Setting up the objective ..." << endl;
    GRBLinExpr exprObj=0;
    for(int b=0; b<nB; b++){
      for(int t=1; t<t_max; t++){
	exprObj += profit[b][t]*(Xvar[b][t]-Xvar[b][t-1]);
      }
      exprObj += profit[b][0]*(Xvar[b][0]);
    }
    //set the objective
    model.setObjective(exprObj,GRB_MAXIMIZE);
    model.update();
    // optimize
    model.optimize();
    // display the solution
    if(show) cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
    // update the upper bound
    upper_bound = model.get(GRB_DoubleAttr_ObjVal);
    // add the constant
    double constant = 0.0;
    for(int r=0; r<r_max; r++){
      for(int t=0; t<t_max; t++){
	char cType = data->getResConstrType(r, t);
	if(cType == 'L'){ 
	  constant += data->getLimit(r, t) * lambda[r][t];
	}
	else if(cType == 'R'){ // not really used for our datasetd
	  constant += lambda[r][t] * (data->getLimit(r, t)); 
	}
	else if(cType == 'I'){// not really used for our datasetd
	  constant -= data->getLimit(r, t, 1);
	  constant += lambda[r][t] * (data->getLimit(r, t)); 
	}
      }
    }
    upper_bound += constant;

    
    // Capture the Y variables    
    for(int b=0; b<nB; b++){
      for(int t = 1; t<t_max; t++){
	if(Xvar[b][t].get(GRB_DoubleAttr_X)-Xvar[b][t-1].get(GRB_DoubleAttr_X)>0.9){
	  int dest = max_block_profit[b][t];
	  Yvar_sol[b][t][dest] = 1;
	}
      }
      if(Xvar[b][0].get(GRB_DoubleAttr_X)>0.9){
	int dest = max_block_profit[b][0];
	Yvar_sol[b][0][dest] = 1;
      }
    }
  } 
  catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } 
  catch(...) {
    cout << "Exception during optimization" << endl;
  }
  if(env) 
    delete env;
  return Yvar_sol;
}

/*
  Functions related to the VLNS
*/

/*
  Reset the lower bounds which were set in the previous iteration
*/

void MIP::resetY(vector<vector<vector<GRBVar> > > &Yvar){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  for(int b=0; b<nB; b++){
    for(int t=0; t<t_max; t++){
      for(int d=0; d<d_max; d++){
	Yvar[b][t][d].set(GRB_DoubleAttr_LB, 0);
	Yvar[b][t][d].set(GRB_DoubleAttr_UB, 1);
      }
    }
  }
}

/*
  Check for resource issues
*/

bool MIP::isResourceFeasible(const vector<vector<vector<double> > > &Yvar_sol){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const int r_max = data->getnResources();
  std::vector<Block> * blocks=data->getBlock();

  bool feas = true;
  for(int r=0; r<r_max; r++){
    for(int t=0; t<t_max; t++){
      long double ru = 0.0;
      for(int b=0; b<nB; b++){
	for(int d=0; d<d_max; d++){
	  double coef = (*blocks)[b].getRCoef(d,r);
	  ru += coef*Yvar_sol[b][t][d];
	}
      }
      if(ru - data->getLimit(r, t) > 0.000001){
	cout << "Problem at resource " << r << ", time " << t << endl;
	cout << "Available: " << data->getLimit(r, t) << ", used: " << ru << endl;
	feas = false;
      }
    }
  }
  return feas;
}



/*
  Solve a MIP, the extended version
  Fix a few variables first
*/

void MIP::MIP_fixed(GRBModel &model, vector<vector<vector<GRBVar> > > &Yvar){

  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  vector<vector<GRBVar> > Xvar(nB);
  try {
    // Create variables and set them to be binary
    // The Y variables are set according to the solution provided
    for(int b=0; b<nB; b++){
      Xvar[b].resize(t_max);
      Yvar[b].resize(t_max);
      for(int t=0; t<t_max; t++){
	Xvar[b][t] = model.addVar(0,1,0,GRB_BINARY);
	Yvar[b][t].resize(d_max);
	for(int d=0; d<d_max; d++){
	  ostringstream str;
	  str << "Y_" << b << "_" << t << "_" << d << endl;
	  //Yvar[b][t][d] = model.addVar(0,1,0,GRB_CONTINUOUS, str.str());
	  // make sure to send everything to a destination or not
	  Yvar[b][t][d] = model.addVar(0,1,0,GRB_BINARY, str.str()); 
	}
      }
    }
    // Update model to integrate the variables
    model.update();
    // type of constraints
    precedenceConstraintsExt(model,Xvar);
    sumDestCompletionConstraint(model,Xvar,Yvar);
    resourceConstraints(model,Yvar);
    // Set up the objective
    setObjective(model,Yvar);
  } 
  catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } 
  catch(...) {
    cout << "Exception during optimization" << endl;
  }
}

/*
  Do a resolve without setting the constraints, etc every time
*/

double MIP::MIP_fixed_resolve(int mod_tmin, int mod_tmax,
			      bool show, GRBModel &model, vector<vector<vector<GRBVar> > > &Yvar,
			      vector<vector<vector<double> > > &Yvar_sol,
			      const vector<double> &profit
			      ){
  
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  double obj = 0.0;

  try {
    // set an initial solution
    // and also fix some variables
    for(int b=0; b<nB; b++){
      for(int t=0; t<t_max; t++){
	for(int d=0; d<d_max; d++){
	  //if((t >= mod_tmin && t <= mod_tmax) && profit[b] < 0.0) 
	  //if((t >= mod_tmin && t <= mod_tmax) || profit[b] < 0.0) 
	  //  continue;
	  if(profit[b] < 0.0) 
	    continue;
	  if(Yvar_sol[b][t][d] == 1.0){
	    Yvar[b][t][d].set(GRB_DoubleAttr_LB, 1);
	  }
	  else if(Yvar_sol[b][t][d] == 0.0){
	    Yvar[b][t][d].set(GRB_DoubleAttr_UB, 0);
	  }
	  // the constraint is dangerous as it could possibly not be undone
	  //model.addConstr(Yvar[b][t][d] == Yvar_sol[b][t][d]);
	}
      }
    }
    // Update model to integrate the variables
    model.update();
    model.optimize();
    // display the solution
    obj = model.get(GRB_DoubleAttr_ObjVal);
    //cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
    for(int b=0; b<nB; b++){
      for(int t =0; t<t_max; t++){
	for(int d =0; d<d_max; d++){
	  Yvar_sol[b][t][d] = 0.0;
	  if(Yvar[b][t][d].get(GRB_DoubleAttr_X) > 0.0){
	    Yvar_sol[b][t][d] = Yvar[b][t][d].get(GRB_DoubleAttr_X);
	  }
	}
      }
    }
    //cout << "Solution" << endl;
  } 
  catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } 
  catch(...) {
    cout << "Exception during optimization" << endl;
  }
  return obj;
}

/*
  Get the average of the negative valued tasks
*/

double MIP::getNegAvgProfit(const vector<double> &profit){
  double sum = 0;
  double length = 0;
  for(int i=0;i<profit.size();i++){
    if(profit[i] < 0){
      sum+=profit[i];
      length++;
    }
  }
  return 2*sum/(4*length);
}

/*
  Check whether the assigned blocks are in the right place
*/

void MIP::checkBlockAssignments(const vector<vector<vector<double> > > &solution, 
				const vector<vector<vector<double> > > &Yvar_sol,
				const vector<double> &profit
				){
  const int nB = data->getNBlock();
  const int d_max = data->getnDestination();

  int done = 0;
  int undone = 0;
  int removed = 0;
  cout << "Negative blocks at this time point: " << endl;
  
  for(int b=0; b<nB; b++){
    for(int t=0; t<1; t++){
      for(int d=0; d<d_max; d++){
	if(solution[b][t][d] > 0.0 && profit[b] < 0.0){
	  cout << "Block: " << b << ", time: " << t 
	       << ", original: " << profit[b] << endl; 
	}
      }
    }
  }
  cout << "Changed blocks ... " << endl;
  for(int b=0; b<nB; b++){
    int fr = 0;
    for(int t=0; t<1; t++){
      for(int d=0; d<d_max; d++){
	if(Yvar_sol[b][t][d] != solution[b][t][d]){
	  cout << "Block: " << b << ", time: " << t 
	       << ", original: " << solution[b][t][d] 
	       << ", new: " << Yvar_sol[b][t][d] << ", profit: " << profit[b] << endl;
	}
	if(Yvar_sol[b][t][d] > solution[b][t][d]){
	  done++;
	  fr++; 
	}
	else if (Yvar_sol[b][t][d] < solution[b][t][d]){
	  undone++;
	  fr--;
	}
      }
    }
    if(fr == -1) removed++;
  }
  cout << "Done: " << done << ", undone: " << undone << ", removed: " << removed << endl;
}

/*
  Update the solution vectors
*/

void MIP::updateSolution(const vector<vector<vector<double> > > &Yvar_sol, 
			 vector<int> &X,
			 vector<vector<double> > &Y
			 ){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();

  X.clear();
  Y.clear();
  X.resize(nB);
  Y.resize(nB);
  for(int b=0; b<nB; b++){
    Y[b].resize(d_max,0.0);
    for(int t=0; t<t_max; t++){
      for(int d=0; d<d_max; d++){
	if(Yvar_sol[b][t][d] > 0.0){
	  X[b] = t;
	  Y[b][d] = Yvar_sol[b][t][d]; 
	}
      }
    }
  }
}

/*
  Update the master solution in parallel
*/

void MIP::updateParSol(MineProblem *problem, const double &obj,
		       vector<vector<vector<double> > > &Yvar_sol
		       ){
  if(problem != NULL){
    if(problem->current->sol_int.obj > obj){ // yes we do
      cout << "\t\tUpdating solution from VLNS (0)..." << endl;
      vector<int> X;// = problem.current->sol_int.x;
      vector<vector<double> > Y;// = problem.current->sol_int.y;
      double obj1;
      int time1;
      problem->get_sol(time1, obj1, X, Y);
      // Now update Yvar_sol correctly
      const int t_max = data->getNPeriod();
      const int nB = data->getNBlock();
      // Do the conversions
      const int d_max = data->getnDestination();  
      for(int b=0; b<nB; b++){
	for(int t=0; t<t_max; t++){
	  for(int d=0; d<d_max; d++){
	    Yvar_sol[b][t][d] = 0.0; // first reset
	    if(X[b] == t && Y[b][d] > 0.0){
	      Yvar_sol[b][t][d] = Y[b][d];
	    }
	  }
	}
      }
    }
  }
}

/*
  Generate a solution needed by the parallel side of things
  Note the objective is not needed here
 */

void MIP::generateParSolution(const vector<vector<vector<double> > > &Yvar_sol,
			      const double &tobj, 
			      Sol_Int &sol){
  const int t_max = data->getNPeriod();
  const int nB = data->getNBlock();
  // Do the conversions
  const int d_max = data->getnDestination();  
  vector<int> X;
  vector<vector<double> > Y;

  X.clear();
  Y.clear();
  X.resize(nB,t_max);
  Y.resize(nB);
  
  for(int b=0; b<nB; b++){
    Y[b].resize(d_max,0.0);
    for(int t=0; t<t_max; t++){
      for(int d=0; d<d_max; d++){
	if(Yvar_sol[b][t][d] > 0.0){
	  X[b] = t;
	  // In the case of the MIP, Y may not be 1
	  Y[b][d] = Yvar_sol[b][t][d];
	}
      }
    }
  }

  sol.x = X;
  sol.y = Y;
  sol.obj = tobj;
  sol.nT = t_max;
}

/*
  Set some default settings for the model parameters
*/

void MIP::setModelParams(GRBModel &model){
  //model.getEnv().set(GRB_IntParam_OutputFlag, 0);
  model.set(GRB_StringAttr_ModelName, "MIP_OPBS");
  model.getEnv().set(GRB_IntParam_BarIterLimit, 30);
  model.getEnv().set(GRB_DoubleParam_MIPGap, 0.005);
  //model.getEnv().set(GRB_DoubleParam_MIPGap, 0.001);
  model.getEnv().set(GRB_IntParam_Presolve, 1);
  model.getEnv().set(GRB_IntParam_Threads, 4);
  model.getEnv().set(GRB_IntParam_BarOrder, 1);
  model.getEnv().set(GRB_IntParam_Cuts,-1); // aggressive cut generation
  model.getEnv().set(GRB_IntParam_CutPasses,200); // maximum number of cutting plane passes for the root cut generation
  //model.getEnv().set(GRB_IntParam_Seed, 10052012);   // Specify only for debugging
  // Not sure what the below mean.
  model.getEnv().set(GRB_DoubleParam_BarConvTol,0.0001);
  model.getEnv().set(GRB_IntParam_PreSparsify,1);
  model.getEnv().set(GRB_IntParam_PrePasses,2);
}

/*
  Update the master solution
  Note: only for the parallel implementation
*/

void MIP::updateMaster(MineProblem *problem,
		  const double &obj,
		  const vector<vector<vector<double> > > &Yvar_sol
		  ){
  if(obj > problem->get_obj()){
    // update problem solution if we have found something better
    // Generally we have with this method
    Sol_Int *sol = new Sol_Int();
    generateParSolution(Yvar_sol,obj,*sol);
    problem->update_sol(*sol);
    cout << "Updated best solution in master  with objective: " << obj << endl;
    problem->set_updated();
    delete sol;
  }
}

double MIP::checkACOSolution( const vector<vector<vector<double> > > &Yvar_sol){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const int r_max = data->getnResources();
  const double rate = data->getDiscountRate();
  std::vector<Block> * blocks=data->getBlock();
  
  vector<int> X(nB,-1);
  vector<vector<double> > Y;
  SolverFunctions *funcs = new SolverFunctions(*data, 210000);
  vector<int> bd(nB);
  for(int b=0; b<nB; b++){
    for(int t =0; t<t_max; t++){
      for(int d =0; d<d_max; d++){
	if (Yvar_sol[b][t][d] == 1.0){
	  X[b] = t;
	  bd[b] = d;
	}
      }
    }
  }
  /*cout << "Checking resources ... " << endl;
  vector<vector<double> > resources_used(t_max);
  for(int t = 0; t < t_max; t++) {
    resources_used[t].resize(r_max,0.0);
  }

  for(int b = 0; b < nB; b++){
    int t = X[b]; 
    if (t==-1) continue;
    int d = bd[b];
    for(int r = 0; r < r_max; r++){
      double total = (*blocks)[b].getRCoef(d,r);
      resources_used[t][r]+=total;
      if(resources_used[t][r] > data->getLimit(r,t)) {
	cout << "Resource problem for block: " << b  << ", used: " 
	     << resources_used[t][r] << ", available: " << data->getLimit(r,t) << endl;
      }
    } 
    double coef = (*blocks)[b].getProfit(d);
    coef /= pow(1+rate,t); 
    //objective += coef;
    } */
  // Construct a solution from problem or X and Y 
  ACO_Solution* sol1 = new ACO_Solution(data,funcs->profit,
					bd,funcs->succ,
					funcs->negBlocks,
					X, Y, true);
  // check the difference between ACO times and MIP times
  vector<int> times = sol1->getBlockTimes();
  /*for(int b=0; b<nB; b++){
    if(times[b] != X[b]){
      cout << "Mismatch for block: " << b << ", aco time: " 
         << times[b] << ", mip time: " << X[b] << endl;
    }
  }
  exit(0);*/

  double obj = sol1->getObj();
  delete sol1;
  delete funcs;
  return obj;
}

/*
  Solve a MIP that is partially fixed
  Similar o the original function with BranchNode_info instead
    note: from bni, only use time here
*/

double MIP::MIP_fixed_resolve(bool show, GRBModel &model, vector<vector<vector<GRBVar> > > &Yvar,
			      vector<vector<vector<double> > > &Yvar_sol,
			      MineProblem *problem
			      ){
  
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();

  //cout << endl << "Fixing and running the extended model ... " << endl;
  double obj = 0.0;

  try {
    // set an initial solution
    // and also fix some variables
    for(int b=0; b<nB; b++){
      if(problem->current->node.time[b][0] == problem->current->node.time[b][1]){ // block can't be mined, set all values to 0
	for(int t=0; t<t_max; t++){
	  for(int d=0; d<d_max; d++){
	    Yvar[b][t][d].set(GRB_DoubleAttr_UB, 0);
	  }
	}
      }
      else{ // otherwise mine from time[b][0] to time[b][1]-1
	for(int t=0; t<t_max; t++){
	  if(t >= problem->current->node.time[b][0] 
	     && t <= problem->current->node.time[b][1]-1)
	    continue;
	  for(int d=0; d<d_max; d++){
	    Yvar[b][t][d].set(GRB_DoubleAttr_UB, 0);
	  }
	}
      }
    }
    // Update model to integrate the variables
    model.update();
    model.optimize();
    // display the solution
    obj = model.get(GRB_DoubleAttr_ObjVal);
    //cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
    for(int b=0; b<nB; b++){
      for(int t =0; t<t_max; t++){
	for(int d =0; d<d_max; d++){
	  Yvar_sol[b][t][d] = 0.0;
	  if(Yvar[b][t][d].get(GRB_DoubleAttr_X) > 0.0){
	    Yvar_sol[b][t][d] = Yvar[b][t][d].get(GRB_DoubleAttr_X);
	  }
	}
      }
    }
    //cout << "Solution" << endl;
  } 
  catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } 
  catch(...) {
    cout << "Exception during optimization" << endl;
  }
  return obj;
}



/*
  Solve a MIP, the extended version
  The time horizon is determined only by mod_tmin an mod_tmax
  blku stores the set of blocks to be optmised during this period
  X is to provide information about blocks that are to be set
  Y is to provide a starting solution
*/

double MIP::MIP_fixed_partial(const vector<bool> &blku,
			      int mod_tmin, int mod_tmax,
			      vector<vector<vector<double> > > &Yvar_sol,
			      const vector<vector<double> > &X,
			      bool last_window
			      ){

  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const int r_max = data->getnResources();
  const double rate = data->getDiscountRate();
  std::vector<Block> * blocks=data->getBlock();

  GRBEnv *env = new GRBEnv();

  vector<vector<GRBVar> > Xvar(nB);
  vector<vector<vector<GRBVar> > > Yvar(nB);
  cout << endl << "\tSetting up the extended model ... " << endl;
  int time = mod_tmax-mod_tmin+1;
  double obj = 0.0;
  try {
    GRBModel model = GRBModel(*env);
    setModelParams(model);
    // Create variables and set them to be binary
    // Only allow positive values for blocks that exist
    for(int b=0; b<nB; b++){
      Xvar[b].resize(time);
      Yvar[b].resize(time);
      for(int t=0; t<time; t++){
	if(blku[b]){
	  if(X[b][t] == 1.0)
	    Xvar[b][t] = model.addVar(1,1,0,GRB_BINARY);
	  else
	    Xvar[b][t] = model.addVar(0,1,0,GRB_BINARY);
	}
	else
	  Xvar[b][t] = model.addVar(0,0,0,GRB_BINARY);
	Yvar[b][t].resize(d_max);
	for(int d=0; d<d_max; d++){
	  ostringstream str;
	  str << "Y_" << b << "_" << t << "_" << d << endl;
	  if(blku[b]){
	    Yvar[b][t][d] = model.addVar(0,1,0,GRB_CONTINUOUS, str.str());
	  }
	  else
	    Yvar[b][t][d] = model.addVar(0,0,0,GRB_CONTINUOUS, str.str());
	}
      }
    }
    // Update model to integrate the variables
    model.update();
    // To set a starting solution
    for(int b=0; b<nB; b++){
      for(int t=0; t<time; t++){
	for(int d=0; d<d_max; d++){
	    Yvar[b][t][d].set(GRB_DoubleAttr_Start, Yvar_sol[b][t+mod_tmin][d]);
	}
      }
    }
    model.update();

    // type of constraints
    ///Precedence (7) 
    for(int a=0; a<nB; a++){
      if(!blku[a]) // a block that is not used
	continue;
      std::vector<int> * pred = (*blocks)[a].getPreds();
      int n = (*blocks)[a].getNumPred();
      for(int p=0; p<n; p++){
	int b = (*pred)[p];
	if(!blku[b]) // a block that is not used
	  continue;
	for(int t=0; t<time; t++){
	  ostringstream str;
	  str << "X_" << a << "_" << t << "_leq_" 
	      << "X_" << b << "_" << t;
	  model.addConstr(Xvar[a][t] <= Xvar[b][t], str.str());
	}
      }
    }
    
    // SumDest (8) 
    for(int b=0; b<nB; b++){
      if(!blku[b]) // a block that is not used
	continue;
      for(int t=1; t<time; t++){
	GRBLinExpr expr = 0;
	expr = Xvar[b][t]-Xvar[b][t-1];
	GRBLinExpr expr2 = 0;
	for(int d=0; d<d_max; d++)
	  expr2 += Yvar[b][t][d];
	model.addConstr(expr2 == expr, "");
      }
      GRBLinExpr expr2 = 0;
      for(int d=0; d<d_max; d++)
	expr2 += Yvar[b][0][d];
      model.addConstr(expr2 == Xvar[b][0], "");
    }
    model.update();

    // changed model, so once a job is completed it stays that way
    for(int b=0; b<nB; b++){
      if(!blku[b]) // a block that is not used
	continue;
      for(int t=0; t<time-1; t++){
	model.addConstr(Xvar[b][t] <= Xvar[b][t+1], "");
      }
    }
    model.update();

    // a new constraint, only needed for this model
    // ensure every block completes
    for(int b=0; b<nB && !last_window; b++){
      if(blku[b])
	model.addConstr(Xvar[b][time-1] == 1, "");
    }
    model.update();
    
    // Resource constraints (10) 
    for(int r=0; r<r_max; r++){
      for(int t=0; t<time; t++){
	char cType = data->getResConstrType(r, t);
	//cout << "cType =" << cType << endl;
	GRBLinExpr expr = 0;
	for(int b=0; b<nB; b++){
	  if(!blku[b]) // a block that is not used
	    continue;
	  for(int d=0; d<d_max; d++){
	    double coef = (*blocks)[b].getRCoef(d,r);
	    expr += coef*Yvar[b][t][d];
	  }
	}
	if(cType == 'L'){ 
	  ostringstream str;
	  str << "Res_" << r << "_" << t;
	  model.addConstr(expr <= data->getLimit(r, t), str.str());
	}
	else if(cType == 'R'){
	  model.addConstr(expr >= data->getLimit(r, t), "");
	}
	else if(cType == 'I'){
	  model.addConstr(expr <= data->getLimit(r, t, 1), "");
	  model.addConstr(expr >= data->getLimit(r,t), "");
	}
      }
    }
    model.update();

    // objective function
    //cout << endl << "Setting up the objective ..." << endl;
    GRBLinExpr exprObj=0;
    for(int b=0; b<nB; b++){
      if(!blku[b]) // a block that is not used
      	continue;
      for(int d=0; d<d_max; d++){
	double b_coef = (*blocks)[b].getProfit(d);
	for(int t=mod_tmin; t<mod_tmax+1; t++){
	  double coef = b_coef/pow(1+rate,t); 
	  exprObj += coef*Yvar[b][t-mod_tmin][d];
	}
      }
    }
    //set the objective and optimise
    model.setObjective(exprObj,GRB_MAXIMIZE);
    model.update();

    model.optimize();
    // display some solution information
    obj = model.get(GRB_DoubleAttr_ObjVal);
    //cout << "\tObj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
    for(int b=0; b<nB; b++){
      //if (!blku[b])
      //continue;
      double b_sum = 0.0;
      for(int t =0; t<time; t++){
	for(int d =0; d<d_max; d++){
	  Yvar_sol[b][t+mod_tmin][d] = Yvar[b][t][d].get(GRB_DoubleAttr_X);
	}
      }
      for(int t =0; t<t_max; t++){
	for(int d =0; d<d_max; d++){
	  b_sum+=Yvar_sol[b][t][d];
	}
      }
      if(b_sum - 1.0 > 0.00001){
	cout << setprecision(20) << "\n\t\tWarning 1:: Block: " << b 
	     << " is used too often, sum: " << b_sum << ", violaton at: " << endl;
	for(int t =0; t<t_max; t++){
	  for(int d =0; d<d_max; d++){
	    if(Yvar_sol[b][t][d] > 0.0){
	      cout << "\t\t\tTime: " << t << ", dest: " << d << endl;
	    }
	  }
	}
      }
    }
  } 
  catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } 
  catch(...) {
    cout << "Exception during optimization" << endl;
  }
  delete env;
  return obj;
}

/*
  A VLNS shell
  Mainly needed to pass a model in without re-creating it each time
  The above will need to be done later.
*/
double MIP::VLNS(const vector<vector<vector<double> > > &solution, 
		 const vector<double> &profit, Timer &timer,
		 vector<int> &X,
		 vector<vector<double> > &Y,
		 MineProblem *problem
		 ){
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  // set up the model now
  GRBEnv *env = new GRBEnv();
  GRBModel model = GRBModel(*env);
  setModelParams(model);
  vector<vector<vector<GRBVar> > > Yvar(nB);
  MIP_fixed(model, Yvar);
  vector<vector<vector<double> > > Yvar_sol = solution;
  int time = 0;
  int shift = 1;
  double obj = 0.0;
  //double treshold = 0.0;//getNegAvgProfit(profit);

  while( time < t_max){ // solve for a time window
    cout << "\n\tFixing blocks before " << time << " and after " << time+shift-1 << endl;
    bool feas = isResourceFeasible(Yvar_sol);
    if(!feas)
      cout << "\t\tWarning: starting with an infeasible solution" << endl;
    // call the MIP
    double t_obj = MIP_fixed_resolve(time, time+shift,true, model, Yvar, Yvar_sol, profit);
    if(t_obj > obj){
      obj = t_obj;
    }
    if(problem != NULL) 
      updateMaster(problem, obj, Yvar_sol);
    // check to see if we have a new best solution from any other thread
    updateParSol(problem,obj,Yvar_sol);
    // reset the model and the Y variables
    model.reset();
    resetY(Yvar);
    cout << "\tCompleted time point " << time << ", CPUT: " 
       << timer.elapsed_time(Timer::VIRTUAL) << ", RT: "  
       << timer.elapsed_time(Timer::REAL)
       << endl;

    time++;
  }  
  // Update the solution vectors
  updateSolution(Yvar_sol,X,Y);
  delete env;
  return obj;
}


/*
  A VLNS shell
  Mainly needed to pass a model in without re-creating it each time
  The above will need to be done later.
*/
double MIP::VLNS_Small(const vector<vector<vector<double> > > &solution, 
		       const vector<double> &profit, Timer &timer,
		       vector<int> &X,
		       vector<vector<double> > &Y,
		       MineProblem *problem
		       ){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const double rate = data->getDiscountRate();
  std::vector<Block> * blocks=data->getBlock();

  vector<vector<vector<double> > > Yvar_sol = solution;
  int time = 0;
  int shift = 1;
  double obj = 0.0;

  //double treshold = 0.0;//getNegAvgProfit(profit);

  // If it is the last wondow, allow (negative) jobs not to complete
  bool last_window = false;

  while( time < t_max-shift){ // solve for a time window
    cout << "\n\tFixing blocks before " << time << " and after " << time+shift-1 << endl;
    bool feas = isResourceFeasible(Yvar_sol);
    if(!feas)
      cout << "\t\tWarning: starting with an infeasible solution" << endl;
    
    double t_obj = 0.0;
    if(time == 0){ // Set this correctly if you need to solve an iteration of the whole lot first
      // set up the model now
      GRBEnv *env = new GRBEnv();
      GRBModel model = GRBModel(*env);
      setModelParams(model);
      vector<vector<vector<GRBVar> > > Yvar(nB);
      MIP_fixed(model, Yvar);
      // call the MIP
      t_obj = MIP_fixed_resolve(time, time+shift, true, model, Yvar, Yvar_sol, profit);
      // Check the difference in the blocks
      //checkBlockAssignments(solution,Yvar_sol,profit);
      //exit(0);
      // Check the ACO solution
      //double obj1 = checkACOSolution(Yvar_sol);
      //cout << "***ACO solution obj: " << obj1 << endl;
      delete env;
    }
    else{
      shift = 2;
      // determine which blocks to be freed during the time points
      vector<bool> blku(nB, false);
      vector<vector<double> > X(nB);
      int b_count = 0;
      for(int b=0; b<nB; b++){
	X[b].resize(shift+1,0.0);
	for(int t=time; t<time+shift+1; t++){
	  for(int d=0; d<d_max; d++){
	    if(Yvar_sol[b][t][d] > 0.0){
	      blku[b] = true;
	      b_count++;
	    }
	    //if(Yvar_sol[b][t][d] > 0.0 && profit[b] > 0.0){
	    if(Yvar_sol[b][t][d] > 0.0){
	      X[b][t-time] = 1.0;
	    }
	  }
	}
      }
      
      cout << "\tFreeing " << b_count << " blocks ..." << endl;
      // call the MIP
      // Note: Yvar_sol will be modified after this call
      t_obj = MIP_fixed_partial(blku, time, time+shift, Yvar_sol, X, last_window);
      // Add the objective values for all blocks at the unfixed time points
      t_obj = 0.0;
      for(int b=0; b<nB; b++){
	double b_sum = 0.0;
	for(int t=0; t<t_max; t++){
	  double sum = 0.0;
	  for(int d=0; d<d_max; d++){
	    if(Yvar_sol[b][t][d] > 0.0){
	      sum += Yvar_sol[b][t][d];
	      b_sum += Yvar_sol[b][t][d];
	      double val = min(Yvar_sol[b][t][d],1.0);
	      double coef = (*blocks)[b].getProfit(d) * val;
	      coef /= pow(1+rate,t); 
	      t_obj += coef;
	    }
	  }
	  //if(sum > 1.0001)
	  //  cout << "\n\t\tWarning:: Block: " << b << ", time: " << t << ", sum: " << sum << endl;
	}
	//if(b_sum > 1.0001)
	//cout << "\n\t\tWarning:: Block: " << b << "is used too often, sum: " << b_sum << endl;
      }
      // check to see if we have a new best solution from any other thread
      updateParSol(problem,obj,Yvar_sol);
      // Check the ACO solution
      //double obj1 = checkACOSolution(Yvar_sol);
      //cout << "***ACO solution obj: " << obj1 << endl;

    }
    cout << "\tRepaired objective: " << t_obj << endl;
    if(t_obj > obj){
      obj = t_obj;
    }
    if(problem != NULL) 
      updateMaster(problem, obj, Yvar_sol);
    feas = isResourceFeasible(Yvar_sol);
    if(!feas){
      cout << "\t\tWarning: solution maybe infeasible, check the resource usage." << endl;
    }
    cout << "\tCompleted time point " << time << ", objective: " << obj << ", CPUT: " 
       << timer.elapsed_time(Timer::VIRTUAL) << ", RT: "  
       << timer.elapsed_time(Timer::REAL)
       << endl;
    time++;
  }
  // Update the solution vectors
  updateSolution(Yvar_sol,X,Y);
  return obj;
}

/*
  A VLNS shell that only work in the window defined by time and time+shift-1
*/
double MIP::VLNS_Window(const vector<vector<vector<double> > > &solution, 
			const vector<double> &profit, Timer &timer,
			vector<int> &X,
			vector<vector<double> > &Y,
			MineProblem *problem,
			int time, int shift
			){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();
  const double rate = data->getDiscountRate();
  std::vector<Block> * blocks=data->getBlock();

  vector<vector<vector<double> > > Yvar_sol = solution;
  //int time = 0;
  //int shift = 1;
  double obj = 0.0;
  // If it is the last wondow, allow (negative) jobs not to complete
  bool last_window = false;
  cout << "\n\tFixing blocks before " << time << " and after " << time+shift-1 << endl;
  double t_obj = 0.0;
  // determine which blocks to be freed during the time points
  vector<bool> blku(nB, false);
  vector<vector<double> > Xt(nB);
  int b_count = 0;
  for(int b=0; b<nB; b++){
    Xt[b].resize(shift+1,0.0);
    for(int t=time; t<time+shift+1; t++){
      for(int d=0; d<d_max; d++){
	if(Yvar_sol[b][t][d] > 0.0){
	  blku[b] = true;
	  b_count++;
	}
	if(Yvar_sol[b][t][d] > 0.0){
	  Xt[b][t-time] = 1.0; // Xt is only the model over time horizon
	}
      }
    }
  }
  
  cout << "\tFreeing " << b_count << " blocks ..." << endl;
  // call the MIP
  // Note: Yvar_sol will be modified after this call
  t_obj = MIP_fixed_partial(blku, time, time+shift, Yvar_sol, Xt, last_window);
  // Add the objective values for all blocks at the unfixed time points
  t_obj = 0.0;
  for(int b=0; b<nB; b++){
    double b_sum = 0.0;
    for(int t=0; t<t_max; t++){
      double sum = 0.0;
      for(int d=0; d<d_max; d++){
	if(Yvar_sol[b][t][d] > 0.0){
	  sum += Yvar_sol[b][t][d];
	  b_sum += Yvar_sol[b][t][d];
	  double val = min(Yvar_sol[b][t][d],1.0);
	  double coef = (*blocks)[b].getProfit(d) * val;
	  coef /= pow(1+rate,t); 
	  t_obj += coef;
	}
      }
    }
  }
  cout << "\tRepaired objective: " << t_obj << endl;
  if(t_obj > obj){
    obj = t_obj;
  }
  cout << "\tCompleted time point " << time << ", objective: " << obj << ", CPUT: " 
       << timer.elapsed_time(Timer::VIRTUAL) << ", RT: "  
       << timer.elapsed_time(Timer::REAL)
       << endl;
  // Update the solution vectors
  updateSolution(Yvar_sol,X,Y);
  return obj;
}

/*
  A VLNS shell (backward), start at time point t_max-1 down to t
  Mainly needed to pass a model in without re-creating it each time
  The above will need to be done later.
*/
double MIP::VLNS_BW(const vector<vector<vector<double> > > &solution, 
		    const vector<double> &profit, Timer &timer,
		    vector<int> &X,
		    vector<vector<double> > &Y,
		    MineProblem *problem
		    ){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();

  // set up the model now
  GRBEnv *env = new GRBEnv();
  GRBModel model = GRBModel(*env);
  setModelParams(model);
  vector<vector<vector<GRBVar> > > Yvar(nB);
  MIP_fixed(model, Yvar);
  vector<vector<vector<double> > > Yvar_sol = solution;
  int shift = 1;
  int time = t_max-shift-1;
  double obj = 0.0;
  //double treshold = 0.0;//getNegAvgProfit(profit);

  while( time > 0){ // solve for a time window
    cout << "\nFixing blocks before " << time << " and after " << time+shift-1 << endl;
    bool feas = isResourceFeasible(Yvar_sol);
    if(!feas)
      cout << "\t\tWarning: starting with an infeasible solution" << endl;
    // set the variables somehow
    // call the MIP
    double t_obj = MIP_fixed_resolve(time, time+shift,true, model, Yvar, Yvar_sol, profit);
    if(t_obj > obj){
      obj = t_obj;
    }
    // Update problem if we have found something better
    if(problem != NULL) 
      updateMaster(problem, obj, Yvar_sol);
    // check to see if we have a new best solution from any other thread
    // Maybe we don't want to do this, so leaving it commented out for now
    //updateParSol(problem);

    // reset the model
    model.reset();
    resetY(Yvar);
    cout << "\tCompleted time point " << time << ", CPUT: " 
       << timer.elapsed_time(Timer::VIRTUAL) << ", RT: "  
       << timer.elapsed_time(Timer::REAL)
       << endl;
    time--;
  }
  // Update the solution vectors
  updateSolution(Yvar_sol,X,Y);
  delete env;
  return obj;
}


/*
  A VLNS shell which fixes 
  Mainly needed to pass a model in without re-creating it each time
  The above will need to be done later.
*/
double MIP::VLNS_fixed(const vector<double> &profit, Timer &timer,
		       vector<int> &X,
		       vector<vector<double> > &Y,
		       MineProblem *problem
		       ){
  // set up some parameters
  const int nB = data->getNBlock();
  const int t_max = data->getNPeriod();
  const int d_max = data->getnDestination();

  // set up the model now
  GRBEnv *env = new GRBEnv();
  GRBModel model = GRBModel(*env);
  setModelParams(model);
  vector<vector<vector<GRBVar> > > Yvar(nB);
  // this can stay the same
  MIP_fixed(model, Yvar);


  // All 0s to being with, should be feasible
  vector<vector<vector<double> > > Yvar_sol(nB);
  for(int b=0; b<nB; b++){
    Yvar_sol[b].resize(t_max);
    for(int t=0; t<t_max; t++){
      Yvar_sol[b][t].resize(d_max,0.0);
    }
  }
  // Rule out some values of Yvar from branch node info
  BranchNode_info bni = problem->get_BranchNode();
  for(int b=0; b<nB; b++){
    int st = bni.time[b][0];
    int et = bni.time[b][1];

    for(int t=0; t<t_max; t++){
      for(int d=0; d<d_max; d++){
	if(st==et){ // can't be mined
	  Yvar[b][t][d].set(GRB_DoubleAttr_UB, 0);
	}
	else if(t < st || t >= et){
	  Yvar[b][t][d].set(GRB_DoubleAttr_UB, 0);
	}
      }
    }
  }
  double obj = 0.0;
  bool show  = true;

  // Nothing to loop over in this variant, maybe for the future?
  //bool feas = isResourceFeasible(Yvar_sol);
  //if(!feas)
  //  cout << "\t\tWarning: starting with an infeasible solution" << endl;
  // call the MIP
  double t_obj = MIP_fixed_resolve(show, model, Yvar, Yvar_sol, problem);
  if(t_obj > obj){
    obj = t_obj;
  }
  // Update problem if we have found something better
  if(problem != NULL) 
    updateMaster(problem, obj, Yvar_sol);
  cout << "\tCompleted run" << ", CPUT: " 
       << timer.elapsed_time(Timer::VIRTUAL) << ", RT: "  
       << timer.elapsed_time(Timer::REAL)
       << endl;

  delete env;
  return obj;
}
// destructor 
MIP::~MIP(){
} 

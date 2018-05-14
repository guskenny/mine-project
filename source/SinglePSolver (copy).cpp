#include "SinglePSolver.h"
#include "SinglePModel.h"

// this is the default one
SinglePSolver::SinglePSolver(int argc,const char **argv){

//  getSettings();
  sh = SettingsHandler();
  sh.printSettings();

  solverType = CPLEX_T; // default is cplex
  solveRelaxed= false;

  // determine solver type and if relaxed
  for(int c=1;c<argc-1;++c){
    std::string flag=argv[c];
    if(flag=="-c")
      solverType = CPLEX_T;
    else if(flag=="-r")
	    solveRelaxed = true;
  }

  // set up new problem model from data file
  std::cout << "setting up model.." << std::flush;
  probModel = new SinglePModel(argv[argc-1]);
  std::cout << " done!" << std::endl;

  std::cout << "setting up max closure factory.. " << std::flush;
  MaxClosureFactory mcfactory;

  for(int i=0;i<argc;++i){
  	if(argv[i][0] == '-' && argv[i][1] != '-')
  	    switch(argv[i][1]){
  		case 'n': case 'b': case 'e': case 'p':
  		case 'N': case 'B': case 'E': case 'P':
  		    mcfactory.setOption(argv[i][1]);
  		    break;
      default:
          mcfactory.setOption('B');
    }
  }

  std::cout << "done!" << std::endl;

  std::cout << "fixing non-upit blocks.. " << std::flush;
  MCSolver = mcfactory(probModel->graph);
  std::cout << "done!" << std::endl << "setting profits.. " << std::flush;
  MCSolver->setProfit(probModel->getProfit());
  std::cout << "done!" << std::endl;
}

SinglePSolver::SinglePSolver(const Daten &prob){

  solverType = CPLEX_T; // default is cplex
  solveRelaxed = false;

  //numPeriods = 2;
  // solve UPIT problem and fix variables
  UpitSolver uSolve = UpitSolver(prob);
  uSolve.solve();
  std::map<int,int> fixed;
  uSolve.getFixed(fixed);

  // create new problem model
  probModel = new SinglePModel(prob);

  MaxClosureFactory mcfactory;
  mcfactory.setOption('B');
  MCSolver = mcfactory(probModel->graph, fixed);
  MCSolver->setProfit(probModel->getProfit());
}

std::vector<double> SinglePSolver::propagateResidual(){//Preprocess &pp){

  std::vector<Block> * blocks=probModel->getBlock();
  const int nB = probModel->getNBlock();
  const int d_max = probModel->getnDestination();

  std::vector<int> succ_count(nB, 0);
  std::vector<double> value(nB, 0.0);

  // initialise succ_count and value vectors
  for (size_t a = 0; a < nB; ++a){
    int numPreds = (*blocks)[a].getNumPred();
    std::vector<int> * pred = (*blocks)[a].getPreds();
    value[a] = residualProfit[a];
    for (size_t p=0; p < numPreds; ++p){
      int b = (*pred)[p];
      succ_count[b]++;
    }
  }

  // find list of blocks with no successors
  std::vector<int> backPropQ;
  for (size_t b = 0; b < nB; ++b){
    if (succ_count[b] == 0)
      backPropQ.push_back(b);
  }

    std::random_device rd; //seed generator
    std::mt19937_64 generator{rd()}; //generator initialized with seed from rd
    //the range is inclusive, so this produces numbers in range [0, 10), same as before
    std::uniform_int_distribution<> dist{0, 99};

  // iterate over backPropQ and randomly distribute among predecessors
  for (size_t a = 0; a < backPropQ.size(); ++a){
    int numPreds = (*blocks)[backPropQ[a]].getNumPred();
    std::vector<int> *pred = (*blocks)[backPropQ[a]].getPreds();

    // generate vector of random numbers to split profits
    std::vector<int> randVec;
    randVec.resize(numPreds);
    double r_sum = 0.0;
    for (size_t i=0; i < randVec.size(); ++i){
      randVec[i] = dist(generator);
      r_sum += randVec[i];
    }
    // scale random vector so it adds up to 1
    std::vector<double> scaledVec;
    scaledVec.resize(randVec.size());
    for (size_t i=0; i < scaledVec.size(); ++i){
      scaledVec[i] = randVec[i]/r_sum;
    }

    for (size_t p = 0; p < numPreds; ++p){
      int b = (*pred)[p];
      value[b] += value[backPropQ[a]] * scaledVec[p];
      if ((--succ_count[b]) == 0)
        backPropQ.push_back(b);
    }
  }

  return value;
}


std::vector<double> SinglePSolver::getFutureProfit(){

  std::vector<Block> * blocks=probModel->getBlock();
  const int nB = probModel->getNBlock();
  const int d_max = probModel->getnDestination();

  std::vector<int> succ_count(nB, 0);
  std::vector<double> value(nB, 0.0);

  // initialise succ_count and value vectors
  for (size_t a = 0; a < nB; ++a){
    int numPreds = (*blocks)[a].getNumPred();
    std::vector<int> * pred = (*blocks)[a].getPreds();
    double maxProfit = 0;
    for (size_t d=0; d < d_max; ++d){
      double profit = (*blocks)[a].getProfit(d);
      if (profit > maxProfit)
        maxProfit = profit;
    }
    value[a] = maxProfit;
    for (size_t p=0; p < numPreds; ++p){
      int b = (*pred)[p];
      succ_count[b]++;
    }
  }

  // find list of blocks with no successors
  std::vector<int> backPropQ;
  for (size_t b = 0; b < nB; ++b){
    if (succ_count[b] == 0)
      backPropQ.push_back(b);
  }

  // iterate over backPropQ
  for (size_t a = 0; a < backPropQ.size(); ++a){
    int numPreds = (*blocks)[backPropQ[a]].getNumPred();
    std::vector<int> *pred = (*blocks)[backPropQ[a]].getPreds();
    double splitProfit = value[backPropQ[a]] / numPreds;
    for (size_t p = 0; p < numPreds; ++p){
      int b = (*pred)[p];
      value[b] += splitProfit;
      if ((--succ_count[b]) == 0)
        backPropQ.push_back(b);
    }
  }

  return value;
}

void SinglePSolver::initModel(qol::MIPSolver &mip,
    std::vector<qol::Variable> &x,
    std::vector<std::vector<qol::Variable> > &y,
    std::vector<Block> * blocks,
    std::vector<double> &profitModifier){

  const int nB = probModel->getNBlock();
  const int r_max = probModel->getnResources();
  const int d_max = probModel->getnDestination();

  std::cout << "setting up decision variables.. " << std::flush;
  // decision variables
  for(int b=0; b<nB; b++){
    x[b]=mip.addVar(0,0,0.0,qol::Variable::BINARY, // LB,UB,cost
        boost::str(boost::format("x%03d")%b));

    y[b].resize(d_max);

    for(int d=0;d<d_max;++d){
      // add modified profit to each destination
      double profit=(*blocks)[b].getProfit(d) + profitModifier[b] * sh.MODIFIER_WEIGHT;
      y[b][d]=mip.addVar(0,0,-profit,qol::Variable::CONTINUOUS,
			      boost::str(boost::format("y%d_%d")%b%d));
      // could write objective as sum expression but this is faster and easy
    }
  }

  std::cout << "done!" << std::endl;
  std::cout << "setting up constraints:" << std::endl << "\tprecedence.. " << std::flush;
  // type of constraints
  /*Precedence (7) (but using different definition of x */
  for(int a=0; a<nB; a++){
    std::vector<int> * pred = (*blocks)[a].getPreds();
    int n = (*blocks)[a].getNumPred();
    for(int p=0; p<n; p++){
      int b = (*pred)[p];
      mip.addConstraint(x[b] >= x[a]).setName("pre%d_%d",b,a);
    }
  }

  std::cout << "done!" << std::endl << "\tsumdest.. " << std::flush;
  /* SumDest (8) */
  for(int b=0; b<nB; b++){
    qol::Expression sumY;
    for(int d=0; d<d_max; d++)
      sumY += y[b][d];
    qol::Expression rhs=x[b];
    mip.addConstraint( sumY == rhs ).setName("dest%d",b);
  }
  std::cout << "done!" << std::endl << "\tresource.. " << std::flush;

  /* Resource constraints (10) */
  for(int r=0; r<r_max; r++){
    char cType = probModel->getResConstrType(r, 0);
    qol::Expression expr;
    for(int b=0; b<nB; b++){
      for(int d=0; d<d_max; d++){
        double coef = (*blocks)[b].getRCoef(d,r);
        if(fabs(coef) > 1e-5)
          expr += coef*y[b][d];
      }
    }
    if(cType == 'L'){
      mip.addConstraint(expr <= probModel->getLimit(r, 0)).setName("R%d",r);
    }else if(cType == 'G'){
      mip.addConstraint(expr >= probModel->getLimit(r, 0)).setName("R%d",r);
    }else{
      std::cerr << "ERROR: resource constraint type " << cType
        << " not implemented - IGNORED\n";
    }
  }
  std::cout << "done!" << std::endl;
}

void SinglePSolver::updateVars(qol::MIPSolver &mip,
  const std::vector<qol::Variable> &x,
  const std::vector<std::vector<qol::Variable> > &y,
  const std::vector<int> &fixed,
  const std::vector<int> &old_fixed,
  const std::vector<std::vector<double> > &y_fixed, int period,
  const BranchNode_info &probInfo,
  const std::vector<bool> &included){

  const double rate = probModel->getDiscountRate();
  const int d_max = probModel->getnDestination();

  int unfixed = 0;

  for (size_t b=0; b < fixed.size(); ++b){

    // unfix all variables in current window that are to be included
    if (included[b] && period >= probInfo.time[b][0] && period < (probInfo.time[b][1]-1)){
      unfixed++;
      if (mip.getVarUB(x[b]) == 0){
        mip.setVarUB(x[b],1);
      }
      for (size_t d=0;d<d_max;++d){
        if (mip.getVarUB(y[b][d]) == 0){
          mip.setVarUB(y[b][d],1);
        }
      }
    }

    if (fixed[b] >= 0 && old_fixed[b] < 0){ // only fix vars that havent been fixed
      mip.setVarLB(x[b],1);
      for (size_t d=0; d < d_max; ++d){
        mip.setVarLB(y[b][d], y_fixed[b][d]);
        //mip.setVarUB(y[b][d], y_fixed[b][d]);
      }
    }
    else{
      if (fixed[b] < 0 && period > 0){ // update profits with discount rate
        for(int d=0;d<d_max;++d){
          double discounted_profit = mip.getObjCoeff(y[b][d]) / (1+rate);
          mip.setObjCoeff(y[b][d], discounted_profit);
        }
      }
    }
  }
  std::cout<< unfixed << " blocks unfixed\n" << std::endl;
}


std::vector<double> SinglePSolver::getConstrUse(const qol::MIPSolver &mip){

  int n_constr = mip.nConstr();
  const int r_max = probModel->getnResources();
  int r = r_max;
  std::vector<double> constrUse (r_max, 0.0);

  for (size_t i=n_constr-1; i > 0; --i){
    qol::ConstraintMIP constr = mip.getConstr(i);
      double prevRHS = mip.getRHS(constr);
      constrUse[--r] = mip.getPrimal(constr);
      if (r == 0)
        break;
  }
  return constrUse;
}


void SinglePSolver::updateConstraints(qol::MIPSolver &mip,
    const std::vector<double> &constrUse,
    const std::vector<double> &prevUse){

  int n_constr = mip.nConstr();
  //std::cout << "n_constr: " << n_constr << std::endl;
  const int r_max = probModel->getnResources();
  int r_count = 0;
  for (size_t i=n_constr-1; i > 0; --i){
    qol::ConstraintMIP constr = mip.getConstr(i);
    std::string name = mip.getConstrName(constr);
    if (name.front() == 'R'){
      //std::cout << "constr" << i << ".name: " << name << std::endl;
      name.erase(name.begin());
      int r = atoi(name.c_str());
      double prevRHS = mip.getRHS(constr);
      //std::cout << "constr" << i << ".RHS = " << prevRHS << std::endl;
      mip.setRHS(constr, prevRHS + (constrUse[r] - prevUse[r]));
      //std::cout << "previous time period use: " << constrUse[r] - prevUse[r] << std::endl;
      //std::cout << "new constr" << i << ".RHS = "
//          << mip.getRHS(constr) << std::endl << std::endl;
      if (r_count++ == r_max)
        break;
    }
  }
}

// function to calculate true objective, as opposed to modified profit
double SinglePSolver::getTrueObj(const Sol_Int &sol){
  const int nB = probModel->getNBlock();
  const int d_max = probModel->getnDestination();
  const double rate = probModel->getDiscountRate();
  std::vector<Block> * blocks=probModel->getBlock();

  double trueObj = 0.0;

  for (size_t b = 0; b < nB; ++b){
    if (sol.x[b] != -1){
      for (size_t d = 0; d < d_max; ++d){
        double blockProfit = (*blocks)[b].getProfit(d);
        trueObj += (blockProfit/pow(1+rate,sol.x[b]))*sol.y[b][d];
      }
    }
  }
  return trueObj;
}

void SinglePSolver::greedySearch(const BranchNode_info &probInfo){
  const int nB = probModel->getNBlock();
  const int t_max = probModel->getNPeriod();
  const int d_max = probModel->getnDestination();
  const int r_max = probModel->getnResources();
  const double rate = probModel->getDiscountRate();

  std::vector<Block> * blocks=probModel->getBlock();

  int block_cnt = 0;

  double mine_total = 0.0;

  std::vector<bool> open_blocks(nB, false);
  std::vector<bool> mined(nB, false);

  for (size_t b=0;b<nB;++b){
    if ((*blocks)[b].getNumPred() == 0 && probInfo.time[b][0] < t_max)
      open_blocks[b] = true;
  }

  bool stop_mining = false;

  int period = 0;

  std::vector<double> res_limit (r_max, 0.0);
  std::vector<double> res_use (r_max, 0.0);

  for (size_t r = 0;r < r_max; ++r){
    res_limit[r] = probModel->getLimit(r, 0);
  }

  while (!stop_mining){
    int max_idx = -1;
    double max_profit = 0;
    int max_dest = 0;
    for (size_t b=0; b<nB;++b){
      // if block open and earliest time less than total time
      if (open_blocks[b] && probInfo.time[b][0] < t_max){
	      const Block & block=probModel->getBlock(b); // get block
        double max_block_profit = block.getProfit(0) / pow(1+rate,period);
        int block_dest = 0;
        for (int d=1;d<d_max;++d){
          double new_profit = block.getProfit(d)/pow(1+rate,period);
          if (new_profit > max_block_profit){
            max_block_profit = new_profit;
            block_dest = d;
          }
        }
        if (max_block_profit > max_profit){
          max_profit = max_block_profit;
          max_idx = b;
          max_dest = block_dest;
        }
      }
    }
    if (max_idx > -1){
      block_cnt++;
      const Block & max_block = probModel->getBlock(max_idx);
      bool change_period = false;
      for (int r = 0; r < r_max; ++r){
        res_use[r] += max_block.getRCoef(max_dest,r);
        if (res_use[r] > res_limit[r]){
          change_period = true;
          break;
        }
      }
      if (change_period){
        for (int r = 0; r < r_max; ++r){
          res_use[r] = max_block.getRCoef(max_dest,r);
        }
        period++;
        if (period >= t_max){
          break;
        }
        max_profit /= (1+rate);
      }
      mine_total += max_profit;
      open_blocks[max_idx] = false;
      mined[max_idx] = true;

      const Node* curr = probModel->graph.getNode(max_idx);
      for (size_t next_arc = 0; next_arc < curr->getOutDegree(); ++next_arc){
        int next_idx = curr->getOutArc(next_arc)->getTgtID();
        const Node* next = probModel->graph.getNode(next_idx);
        bool block_open = true;
        for (size_t prev_arc = 0; prev_arc < next->getInDegree(); ++prev_arc){
          int prev_idx = next->getInArc(prev_arc)->getSrcID();
          if (open_blocks[prev_idx]){
            block_open = false;
            break;
          }
        }
        if (block_open && !mined[next_idx] && probInfo.time[next_idx][0] < t_max){
          open_blocks[next_idx] = true;
        }
      }
    }
    else {
      std::cout << "stop!" << std::endl;
      stop_mining = true;
    }
  }

  std::cout << "\nGreedy search completed!\n\n" << block_cnt << " blocks mined\n"
            << mine_total << " total mine value\n\n";


 return mine_total;
}

// set up MIP and solve
int SinglePSolver::solve(){

  qol::CpuTimer timer;

  qol::MIPSolver *mipPtr=0; // pointer to qol MIP solver

  const int nB = probModel->getNBlock();
  const int t_max = probModel->getNPeriod();
  const int d_max = probModel->getnDestination();
  const int r_max = probModel->getnResources();
  std::vector<Block> * blocks=probModel->getBlock();

  int n_periods = t_max;

  std::string lpFile="";// ,solnFile = probModel->getName()+".sol";
  Sol_Int best_sol(nB, d_max);
  double bestObj = -1e9;

  // establish vector for solutions
  std::vector<Sol_Int> sols;

  if (sh.TEST_SOL_MERGE){
    std::cout << "Testing SolutionMerger class!" << std::endl;

    SolutionMerger sm(sh);

    Sol_Int mergeTest;

    mergeTest.init(4,2);
    mergeTest.x = std::vector<int>{2,3,1,0};
    mergeTest.nT = 5;
    sols.push_back(mergeTest);
    mergeTest.init(4,2);
    mergeTest.x = std::vector<int>{2,4,3,2};
    mergeTest.nT = 5;
    sols.push_back(mergeTest);
    mergeTest.init(4,2);
    mergeTest.x = std::vector<int>{0,2,2,1};
    mergeTest.nT = 5;
    sols.push_back(mergeTest);

    BranchNode_info merged = sm.mergePCPSP(sols);

    std::cout << std::endl << "verified solutions:" << std::endl;

    for (size_t i = 0; i<merged.time[0].size();++i){
      for (size_t b=0; b<merged.time.size(); ++b){
        std::cout << merged.time[b][i] << " ";
      }
      std::cout << std::endl;
    }

    std::cout << std::endl << "testing copy:" << std::endl;

    Sol_Int mergeTest2 = mergeTest;

    std::cout << "mergeTest: ";
    for (int i=0;i<mergeTest.x.size();++i){
      std::cout << mergeTest.x[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "mergeTest2: ";
    for (int i=0;i<mergeTest2.x.size();++i){
      std::cout << mergeTest2.x[i] << " ";
    }
    std::cout << std::endl << std::endl << "testing change:" << std::endl;
    mergeTest2.x = std::vector<int>{1,1,1,1};
    std::cout << "mergeTest: ";
    for (int i=0;i<mergeTest.x.size();++i){
      std::cout << mergeTest.x[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "mergeTest2: ";
    for (int i=0;i<mergeTest2.x.size();++i){
      std::cout << mergeTest2.x[i] << " ";
    }
    std::cout << std::endl;

    return 0;
  }

  for (int iter = 0; iter < sh.NUM_ITER; ++iter){
    std::cout << "\n******* ITERATION " << iter << " of "
              << sh.NUM_ITER << " *******\n" << std::endl;
  double trueObj = 0;
  Sol_Int test_sol(nB, d_max);
   try{
    // set verbosity
    qol::Parameters param;
    param.setParamVal(qol::VERBOSITY,1);

    // vector to maintain fixed blocks
    std::vector<int> fixed(nB, -1); // -1 means not fixed
    std::vector<int> old_fixed(nB, -1);
    std::vector<std::vector<double> > y_fixed(nB, std::vector<double> (d_max, -1));

    std::vector<bool> mined (nB,false);

    // determine solver type
    //if(solverType == GUROBI_T)
    //  mipPtr = new qol::GurobiFormulation();
    //else if(solverType == CPLEX_T)
    if(solverType == CPLEX_T)
	    mipPtr = new qol::CplexFormulation();

    // create MIPSolver object
    qol::MIPSolver &mip=*mipPtr;
    mip.setParameters(param);

    std::vector<qol::Variable> x(nB);
    std::vector<std::vector<qol::Variable> > y(nB);

    std::cout << "Problem has " << nB << " blocks" << std::endl;

    BranchNode_info probInfo(nB,t_max,d_max,r_max);

    Preprocess pp((*probModel),probInfo,mined,sh,true);

//    std::vector<bool> processable;
//
//    if (sh.PROCESS_THRESHOLD > 0){
//      int p_count = pp.getProcessable(processable);
//      std::cout << p_count << " blocks identified as unprocessable" << std::endl;
//    }
//    else
//      processable.resize(nB, true);



    // propagate the residual profit to create a profit modifier
    //std::cout << "propagating residual.." << std::flush;
    std::vector<double> profitModifier(nB,0.0);
    //profitModifier = propagateResidual(pp);
    //std::cout << " done!" << std::endl;

    if (sh.GREEDY_SEARCH){
      n_periods = 0;
      greedySearch(probInfo);
    }
    else{
      initModel(mip, x, y, blocks, profitModifier);
    }
    double obj_value = 0;

    // vector to track previous time period use of resources
    std::vector<double> prevUse(r_max, 0.0);


    for (int t = 0; t < n_periods; ++t){

      std::cout << "\n*** TIME POINT " << t << " ***\n" << std::endl;

      qol::CpuTimer sub_timer;
      // update fixEarliest
      //if (t > 0)

      pp.fixEarliest(mined,t);

      std::vector<bool> included;
      pp.getMostValuableCones(t,included,mined);


      std::cout << std::endl << sub_timer.elapsedSeconds() << " sec CPU / "
                << sub_timer.elapsedWallTime() << " sec (wall) spent pre-processing"
                << std::endl << std::endl;

      std::vector<double> constrUse;

      if (t > 0){
        constrUse = getConstrUse(mip);
      }

      updateVars(mip, x, y, fixed, old_fixed, y_fixed, t, probInfo,included);

      old_fixed = fixed; // update old fixed list

      if (t > 0){
        updateConstraints(mip, constrUse, prevUse);
        //updateConstraints(mip);
      }

      prevUse = constrUse;

      if(lpFile != ""){
        mip.writeLP(lpFile.c_str());
        std::cout << "Wrote " << lpFile << std::endl;
      }

      qol::Status status = solveRelaxed ? mip.solveRelaxed() : mip.solveExact();
      std::cout << boost::format("Completed in %.2f sec CPU / %.2f sec wall. Objective = %f\n"
  			       ) % timer.elapsedSeconds() % timer.elapsedWallTime() % (-mip.getObjective());

      obj_value = -mip.getObjective();

      std::string solnFile = "solfiles/" + probModel->getName()+"_"+std::to_string(t)+".sol";

      std::ofstream fp_out;
      fp_out.open(solnFile.c_str()); //, std::ios::app);

      size_t block_count = 0;
      std::vector<double> d_count(d_max,0.0);

      if(fp_out.is_open()){
        std::cout << "Writing solution to " << solnFile << std::endl;
        fp_out<< "# Status          "<<status<< std::endl;
        fp_out<< "# CPU time        "<<timer.elapsedSeconds()
          <<" sec.  Wall:"<< timer.elapsedWallTime() <<std::endl;
        fp_out<< "# Objective Value "<< obj_value << std::endl;

        fp_out<< "# Block destination time yval"<<std::endl;
        /* solution file format blocks -> destinations -> time non-zero y value*/
        for(int b=0; b<y.size(); b++){
          bool already_counted = false;
  	      for(int d=0;d<d_max;++d){
  	        if(mip.getPrimal(y[b][d]) > 1e-5){
              mined[b] = true;
              if (old_fixed[b] < 0){
                d_count[d] += mip.getPrimal(y[b][d]);
                if (!already_counted){
                  block_count++;
                  already_counted=true;
                }
                fixed[b] = t;
                y_fixed[b][d] = mip.getPrimal(y[b][d]);
              }
  	          fp_out << b << " " << d << " " << fixed[b] << " "
                << mip.getPrimal(y[b][d]) << std::endl;
           }
  	      }
        }
        fp_out.close();
      }

      std::cout << "\nMined " << block_count << " blocks in period " << t << std::endl;

      for (size_t d=0; d<d_max;++d)
        std::cout << d_count[d] << " blocks sent to destination " << d << std::endl;


      // test solution for feasibility
      test_sol.init(nB, d_max);
      int max_period = 0;

      for (size_t b = 0; b < fixed.size(); ++b){
        test_sol.x[b] = fixed[b];
        if (fixed[b] > max_period){
          max_period = fixed[b];
        }
        for (size_t d=0; d < d_max; ++d){
    	    if(mip.getPrimal(y[b][d]) > 1e-5){
            test_sol.y[b][d] = mip.getPrimal(y[b][d]);
          }
        }
      }

      trueObj = getTrueObj(test_sol);

      //test_sol.obj = obj_value;
      test_sol.obj = trueObj;

      test_sol.nT = max_period + 1;

      std::cout << "\nTrue objective: " << trueObj << std::endl;

      std::cout << "\nChecking solution from solver...\n" << std::endl;

      bool test_error = verify((*probModel), test_sol);

      std::cout << "Solution found by solver was ";
      if (!test_error)
        std::cout << "feasible!\n" << std::endl;
      else{
        std::cout << "infeasible!\n" << std::endl;
        throw qol::Exception("Infeasible solution!");
      }

    } // end for loop
      if (trueObj > bestObj){
        std::cout << std::endl << "***** NEW BEST OBJECTIVE FOUND! *****\n" << std::endl;
        bestObj = trueObj;
        best_sol = test_sol;
      }


   } // end try statement
  catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
  catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}


   }//end of iteration loop
    std::vector<std::vector<std::vector<double> > > Y_sol;
    std::cout << "Extracting solution from Sol_Int..." << std::endl;
    extractSol(best_sol,Y_sol);

    std::cout << "*** testing pack/extract functions ***\n" << std::endl;

    Sol_Int packedSol(nB, d_max);

    packSol(Y_sol, packedSol);

    double packedObj = getTrueObj(packedSol);

    std::cout << "packed/extracted obj: " << packedObj << std::endl;

    double improvedObj = 0.0;

    if (sh.VNS_WINDOW_LIMIT){
      improvedObj = windowVNS(Y_sol);
    }
    else if (sh.RECORD_PASS_TIME){
      improvedObj = timeWindow(Y_sol);
    }
    else{
      improvedObj = slideWindow(Y_sol);
    }

    std::cout << "\noriginal objective: " << bestObj << std::endl;
    std::cout << "\nimproved objective: " << improvedObj << std::endl;

    std::cout << "# CPU time        "<<timer.elapsedSeconds()
      <<" sec.  Wall:"<< timer.elapsedWallTime() <<std::endl;

    sh.printSettings();

		// save objective for current run to file (for batch processing)
    if (sh.RECORD_RUN){
      std::string outobjfile = "objfiles/" + probModel->getName()+"_outobj.csv";
      std::ofstream outobj;
		  outobj.open(outobjfile, std::ios_base::app);
		  outobj << bestObj << "," << improvedObj << ","
             << timer.elapsedWallTime() << std::endl;
      std::cout << "objectives and times written to " << outobjfile << std::endl << std::endl;
    }
    std::cout << "# CPU time        "<<timer.elapsedSeconds()
          <<" sec.  Wall:"<< timer.elapsedWallTime() <<std::endl;

  delete mipPtr;
  return 0;
}

double SinglePSolver::SinglePSolve(const BranchNode_info &init_sol, Sol_Int &sol){
  qol::CpuTimer timer;

  qol::MIPSolver *mipPtr=0; // pointer to qol MIP solver

  const int nB = probModel->getNBlock();
  const int t_max = probModel->getNPeriod();
  const int d_max = probModel->getnDestination();
  const int r_max = probModel->getnResources();
  std::vector<Block> * blocks=probModel->getBlock();

  int n_periods = t_max;

  std::string lpFile="";// ,solnFile = probModel->getName()+".sol";

  double bestObj = -1e9;

//  for (int iter = 0; iter < sh.NUM_ITER; ++iter){
//    std::cout << "\n******* ITERATION " << iter << " of "
//              << sh.NUM_ITER << " *******\n" << std::endl;
  double trueObj = 0;

   try{
    // set verbosity
    qol::Parameters param;
    param.setParamVal(qol::VERBOSITY,1);

    // vector to maintain fixed blocks
    std::vector<int> fixed(nB, -1); // -1 means not fixed
    std::vector<int> old_fixed(nB, -1);
    std::vector<std::vector<double> > y_fixed(nB, std::vector<double> (d_max, -1));

    std::vector<bool> mined (nB,false);

    // determine solver type
    //if(solverType == GUROBI_T)
    //  mipPtr = new qol::GurobiFormulation();
    //else if(solverType == CPLEX_T)
    if(solverType == CPLEX_T)
	    mipPtr = new qol::CplexFormulation();

    // create MIPSolver object
    qol::MIPSolver &mip=*mipPtr;
    mip.setParameters(param);

    std::vector<qol::Variable> x(nB);
    std::vector<std::vector<qol::Variable> > y(nB);

    std::cout << "Problem has " << nB << " blocks" << std::endl;

    BranchNode_info probInfo = init_sol;

    Preprocess pp((*probModel),probInfo,mined,sh,true);

//    std::vector<bool> processable;
//
//    if (sh.PROCESS_THRESHOLD > 0){
//      int p_count = pp.getProcessable(processable);
//      std::cout << p_count << " blocks identified as unprocessable" << std::endl;
//    }
//    else
//      processable.resize(nB, true);



    // propagate the residual profit to create a profit modifier
    //std::cout << "propagating residual.." << std::flush;
    std::vector<double> profitModifier(nB,0.0);
    //profitModifier = propagateResidual(pp);
    //std::cout << " done!" << std::endl;

    if (sh.GREEDY_SEARCH){
      n_periods = 0;
      greedySearch(probInfo);
    }
    else{
      initModel(mip, x, y, blocks, profitModifier);
    }
    double obj_value = 0;

    // vector to track previous time period use of resources
    std::vector<double> prevUse(r_max, 0.0);


    for (int t = 0; t < n_periods; ++t){

      std::cout << "\n*** TIME POINT " << t << " ***\n" << std::endl;

      qol::CpuTimer sub_timer;
      // update fixEarliest
      //if (t > 0)

      pp.fixEarliest(mined,t);

      std::vector<bool> included;
      pp.getMostValuableCones(t,included,mined);


      std::cout << std::endl << sub_timer.elapsedSeconds() << " sec CPU / "
                << sub_timer.elapsedWallTime() << " sec (wall) spent pre-processing"
                << std::endl << std::endl;

      std::vector<double> constrUse;

      if (t > 0){
        constrUse = getConstrUse(mip);
      }

      updateVars(mip, x, y, fixed, old_fixed, y_fixed, t, probInfo,included);

      old_fixed = fixed; // update old fixed list

      if (t > 0){
        updateConstraints(mip, constrUse, prevUse);
        //updateConstraints(mip);
      }

      prevUse = constrUse;

      if(lpFile != ""){
        mip.writeLP(lpFile.c_str());
        std::cout << "Wrote " << lpFile << std::endl;
      }

      qol::Status status = solveRelaxed ? mip.solveRelaxed() : mip.solveExact();
      std::cout << boost::format("Completed in %.2f sec CPU / %.2f sec wall. Objective = %f\n"
  			       ) % timer.elapsedSeconds() % timer.elapsedWallTime() % (-mip.getObjective());

      obj_value = -mip.getObjective();

      std::string solnFile = "solfiles/" + probModel->getName()+"_"+std::to_string(t)+".sol";

      std::ofstream fp_out;
      fp_out.open(solnFile.c_str()); //, std::ios::app);

      size_t block_count = 0;
      std::vector<double> d_count(d_max,0.0);

      if(fp_out.is_open()){
        std::cout << "Writing solution to " << solnFile << std::endl;
        fp_out<< "# Status          "<<status<< std::endl;
        fp_out<< "# CPU time        "<<timer.elapsedSeconds()
          <<" sec.  Wall:"<< timer.elapsedWallTime() <<std::endl;
        fp_out<< "# Objective Value "<< obj_value << std::endl;

        fp_out<< "# Block destination time yval"<<std::endl;
        /* solution file format blocks -> destinations -> time non-zero y value*/
        for(int b=0; b<y.size(); b++){
          bool already_counted = false;
  	      for(int d=0;d<d_max;++d){
  	        if(mip.getPrimal(y[b][d]) > 1e-5){
              mined[b] = true;
              if (old_fixed[b] < 0){
                d_count[d] += mip.getPrimal(y[b][d]);
                if (!already_counted){
                  block_count++;
                  already_counted=true;
                }
                fixed[b] = t;
                y_fixed[b][d] = mip.getPrimal(y[b][d]);
              }
  	          fp_out << b << " " << d << " " << fixed[b] << " "
                << mip.getPrimal(y[b][d]) << std::endl;
           }
  	      }
        }
        fp_out.close();
      }

      std::cout << "\nMined " << block_count << " blocks in period " << t << std::endl;

      for (size_t d=0; d<d_max;++d)
        std::cout << d_count[d] << " blocks sent to destination " << d << std::endl;


      // test solution for feasibility
      sol.init(nB, d_max);
      int max_period = 0;

      for (size_t b = 0; b < fixed.size(); ++b){
        sol.x[b] = fixed[b];
        if (fixed[b] > max_period){
          max_period = fixed[b];
        }
        for (size_t d=0; d < d_max; ++d){
    	    if(mip.getPrimal(y[b][d]) > 1e-5){
            sol.y[b][d] = mip.getPrimal(y[b][d]);
          }
        }
      }

      trueObj = getTrueObj(sol);

      //test_sol.obj = obj_value;
      sol.obj = trueObj;

      sol.nT = max_period + 1;

      std::cout << "\nTrue objective: " << trueObj << std::endl;

      std::cout << "\nChecking solution from solver...\n" << std::endl;

      bool test_error = verify((*probModel), sol);

      std::cout << "Solution found by solver was ";
      if (!test_error)
        std::cout << "feasible!\n" << std::endl;
      else{
        std::cout << "infeasible!\n" << std::endl;
        throw qol::Exception("Infeasible solution!");
      }

   } // end try statement
  catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
  catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

  delete mipPtr;
  return trueObj;
}

void SinglePSolver::extractSol(const Sol_Int &sol, std::vector<std::vector<std::vector<double> > > &Y_sol){

  const int nB = probModel->getNBlock();
  const int t_max = probModel->getNPeriod();
  const int d_max = probModel->getnDestination();

  // reset Y solution vector
  Y_sol.clear();
  Y_sol.resize(nB,
      std::vector<std::vector<double> >(t_max,
        std::vector<double>(d_max, 0.0)));

  for (size_t b=0; b<nB; ++b){
    if (sol.x[b] >= 0){
      for (size_t d=0; d<d_max; ++d){
        Y_sol[b][sol.x[b]][d] = sol.y[b][d];
      }
    }
  }
}

void SinglePSolver::packSol(const std::vector<std::vector<std::vector<double> > > &Y_sol, Sol_Int &sol){
  const int nB = probModel->getNBlock();
  const int t_max = probModel->getNPeriod();
  const int d_max = probModel->getnDestination();


  sol.init(nB, d_max);

  for (size_t b=0;b<nB;++b){
    bool mined = false;
    for (size_t t=0;t<t_max;++t){
      for (size_t d=0;d<d_max;++d){
        if (Y_sol[b][t][d] > 0.0){
          sol.y[b][d] = Y_sol[b][t][d];
          sol.x[b] = t;
          mined = true;
        }
      }
    }
    if (!mined)
      sol.x[b] = -1;
  }

  sol.nT = t_max;
  sol.obj = getTrueObj(sol);

}

double SinglePSolver::windowVNS(std::vector<std::vector<std::vector<double> > > &Y_sol){

  qol::CpuTimer timer;

  const int t_max = probModel->getNPeriod();
  const double eps = 1e-6;

  bool early_finish = false;

  std::cout << "\n*** running windowVNS() ***\n"<< std::endl;

  std::cout << "Start window size: " << sh.WINDOW_SIZE << std::endl;
  std::vector<double> objectives;
  std::vector<double> times;

  int window_size = sh.WINDOW_SIZE;


  double obj = 0.0;
  double old_obj = -1e9;

	int converge_count = 0;

  for (size_t i=0;i<sh.NUM_PASSES && !early_finish;++i){

    std::cout << "\n*** running backward pass "<< i+1 << " ***" << std::endl;

    for (size_t t=t_max-window_size-1; t>=0 && !early_finish; --t){
      obj = improveWindow(Y_sol, t, i, window_size, false);
      std::cout << "\nImproved objective: " << obj << std::endl;

      if (sh.OBJ_UB > 0)
        std::cout << "\nObjective upper bound: " << sh.OBJ_UB << std::endl;

      objectives.push_back(obj);
      times.push_back(timer.elapsedWallTime());

      if (sh.RECORD_DATA){
        std::ofstream outfile;
        outfile.open ("timeobj.csv");
        for (size_t i=0; i<objectives.size(); ++i)
          outfile << times[i] << "," << objectives[i] << std::endl;
        outfile.close();
      }

      std::cout << std::endl << "WINDOW VNS TOTAL TIME: " << timer.elapsedSeconds()
                << " sec CPU / " << timer.elapsedWallTime() << " sec wall "
                << std::endl;

      if (sh.WALL_SEARCH_TIME > 0)
        std::cout << "\nTime limit (wall): " << sh.WALL_SEARCH_TIME << std::endl;

      if (sh.CPU_SEARCH_TIME > 0)
        std::cout << "\nTime limit (CPU): " << sh.CPU_SEARCH_TIME << std::endl;

      Sol_Int test_sol;
      packSol(Y_sol, test_sol);

      test_sol.nT = t_max;

      std::cout << "\nChecking solution from solver...\n" << std::endl;

      bool test_error = verify((*probModel), test_sol);

      std::cout << "Solution found by solver was ";
      if (!test_error)
        std::cout << "feasible!\n" << std::endl;
      else{
        std::cout << "infeasible!\n" << std::endl;
        throw qol::Exception("Infeasible solution!");
      }

      if (obj > sh.OBJ_UB + eps && sh.OBJ_UB > 0){
        std::cout << "\nobjective: " << obj << " over upper bound of "
                  << sh.OBJ_UB << ", early finish!!" << std::endl;
        early_finish=true;
      }

      if (sh.WALL_SEARCH_TIME > 0 && timer.elapsedWallTime() > sh.WALL_SEARCH_TIME){
        std::cout << std::endl << "wall time limit of " << sh.WALL_SEARCH_TIME
                  << " exceeded, early finish!!" << std::endl;
        early_finish=true;
      }
      if (sh.CPU_SEARCH_TIME > 0 && timer.elapsedSeconds() > sh.CPU_SEARCH_TIME){
        std::cout << std::endl << "CPU time limit of " << sh.CPU_SEARCH_TIME
                  << " exceeded, early finish!!" << std::endl;
        early_finish=true;
      }
    }

    if (obj - old_obj > eps){
      old_obj = obj;
			converge_count = 0;
      if (window_size > sh.WINDOW_SIZE){
        std::cout << "********** IMPROVEMENT FOUND! RETURNING "
                << "WINDOW SIZE TO " << sh.WINDOW_SIZE
                << " **********\n"<< std::endl;
        window_size = sh.WINDOW_SIZE;
      }
    }
    else{
      converge_count++;
        std::cout << "********** NO IMPROVEMENT IN " << converge_count
                << " FULL PASSES **********\n"<< std::endl;
		}

		if (sh.CONVERGE_STOP && converge_count >= sh.CONVERGE_STOP) {
      window_size++;
      converge_count = 0;

      if (window_size > sh.VNS_WINDOW_LIMIT || window_size > t_max){
        std::cout << "********** WINDOW SIZE LIMIT REACHED! ENDING RUN!"
                  << " **********\n"<< std::endl;
        early_finish = true;
      }
      else{
        std::cout << "********** CONVERGE LIMIT REACHED "
                  << " INCREASING WINDOW SIZE TO " << window_size
                  << " **********\n"<< std::endl;
      }
    }

    if (!early_finish)
      std::cout << "*** running foward pass " << i+1 << " ***" << std::endl;

    for (size_t t=0; t<t_max-sh.WINDOW_SIZE && !early_finish; ++t){
      obj = improveWindow(Y_sol, t, i, window_size, true);
      std::cout << "\nImproved objective: " << obj << std::endl;

      if (sh.OBJ_UB > 0)
        std::cout << "\nObjective upper bound: " << sh.OBJ_UB << std::endl;

      objectives.push_back(obj);
      times.push_back(timer.elapsedWallTime());

      if (sh.RECORD_DATA){
        std::ofstream outfile;
        outfile.open ("timeobj.csv");
        for (size_t i=0; i<objectives.size(); ++i)
          outfile << times[i] << "," << objectives[i] << std::endl;
        outfile.close();
      }

      std::cout << std::endl << "WINDOW VNS TOTAL TIME: " << timer.elapsedSeconds()
                << " sec CPU / " << timer.elapsedWallTime() << " sec wall "
                << std::endl;

      if (sh.WALL_SEARCH_TIME > 0)
        std::cout << "\nTime limit (wall): " << sh.WALL_SEARCH_TIME << std::endl;

      if (sh.CPU_SEARCH_TIME > 0)
        std::cout << "\nTime limit (CPU): " << sh.CPU_SEARCH_TIME << std::endl;


      std::cout << "\nChecking solution from solver...\n" << std::endl;

      Sol_Int test_sol;
      packSol(Y_sol, test_sol);

      test_sol.nT = t_max;

      bool test_error = verify((*probModel), test_sol);

      std::cout << "Solution found by solver was ";
      if (!test_error)
        std::cout << "feasible!\n" << std::endl;
      else{
        std::cout << "infeasible!\n" << std::endl;
        throw qol::Exception("Infeasible solution!");
      }

      if (obj > sh.OBJ_UB + eps && sh.OBJ_UB > 0){
        std::cout << "\nobjective: " << obj << " over upper bound of "
                  << sh.OBJ_UB << ", early finish!!" << std::endl;
        early_finish=true;
      }

      if (sh.WALL_SEARCH_TIME > 0 && timer.elapsedWallTime() > sh.WALL_SEARCH_TIME){
        std::cout << std::endl << "wall time limit of " << sh.WALL_SEARCH_TIME
                  << " exceeded, early finish!!" << std::endl;
        early_finish=true;
      }

      if (sh.CPU_SEARCH_TIME > 0 && timer.elapsedSeconds() > sh.CPU_SEARCH_TIME){
        std::cout << std::endl << "CPU time limit of " << sh.CPU_SEARCH_TIME
                  << " exceeded, early finish!!" << std::endl;
        early_finish=true;
      }
    }

    if (!early_finish){
      if (obj - old_obj > eps){
        old_obj = obj;
  			converge_count = 0;
        if (window_size > sh.WINDOW_SIZE){
          std::cout << "********** IMPROVEMENT FOUND! RETURNING "
                    << "WINDOW SIZE TO " << sh.WINDOW_SIZE
                    << " **********\n"<< std::endl;
          window_size = sh.WINDOW_SIZE;
        }
      }
      else{
        converge_count++;
        std::cout << "********** NO IMPROVEMENT IN " << converge_count
                  << " FULL PASSES **********\n"<< std::endl;

  		}

  		if (sh.CONVERGE_STOP && converge_count >= sh.CONVERGE_STOP) {
        window_size++;
        converge_count=0;

        if (window_size > sh.VNS_WINDOW_LIMIT || window_size > t_max){
          std::cout << "********** WINDOW SIZE LIMIT REACHED! ENDING RUN!"
                    << " **********\n"<< std::endl;
          early_finish = true;
        }
        else{
          std::cout << "********** CONVERGE LIMIT REACHED "
                    << " INCREASING WINDOW SIZE TO " << window_size
                    << " **********\n"<< std::endl;
        }

      }
    }
  }

  if (sh.RECORD_DATA)
    std::cout << "\nslideWindow() complete! Data available in timeobj.csv" << std::endl;
  else{
    std::cout << "\nslideWindow() complete. time, objective values:" << std::endl;
    for (size_t i=0; i<objectives.size(); ++i)
      std::cout << times[i] << ", " << objectives[i] << std::endl;
  }
  std::cout << std::endl;

  return obj;
}

double SinglePSolver::timeWindow(std::vector<std::vector<std::vector<double> > > &Y_sol){

  const int t_max = probModel->getNPeriod();
  const int d_max = probModel->getnDestination();
  const int nB = probModel->getNBlock();
  const double eps = 1e-6;

  bool early_finish = false;

  std::cout << "\n*** running timeWindow() ***\n"<< std::endl;

  std::cout << "Window size: " << sh.WINDOW_SIZE << std::endl;
  std::vector<double> objectives;
  std::vector<double> times;

  Sol_Int orig_sol;
  packSol(Y_sol, orig_sol);

  double orig_obj = getTrueObj(orig_sol);

  double obj = 0.0;
  double old_obj = -1e9;

	int converge_count = 0;

  bool dir_indicator = false;

  for (size_t tests = 0; tests < 2; ++tests){

  for (size_t i=0;i<sh.NUM_PASSES ;++i){

    std::cout << "\n*** running backward pass "<< i+1
      << " WS: " << sh.WINDOW_SIZE << " ***" << std::endl;

    std::vector<std::vector<std::vector<double> > >
        tempY_sol(nB,std::vector<std::vector<double > >(t_max,
              std::vector<double> (d_max,0.0)));


    for (int yb = 0; yb < nB; ++yb){
      for (int yt = 0; yt < t_max; ++yt){
        for (int yd = 0; yd < d_max; ++yd){
          tempY_sol[yb][yt][yd] = Y_sol[yb][yt][yd];
        }
      }
    }

    qol::CpuTimer timer;

    for (int t=t_max-sh.WINDOW_SIZE-1; t>=0; --t){

      std::cout << "** CURRENT TIME: " << t << "**" << std::endl;
      obj = improveWindow(tempY_sol, t, i, sh.WINDOW_SIZE, dir_indicator);
      std::cout << "\nImproved objective: " << obj << std::endl;

      Sol_Int test_sol;
      packSol(Y_sol, test_sol);

      test_sol.nT = t_max;

      std::cout << "\nChecking solution from solver...\n" << std::endl;

      bool test_error = verify((*probModel), test_sol);

      std::cout << "Solution found by solver was ";
      if (!test_error)
        std::cout << "feasible!\n" << std::endl;
      else{
        std::cout << "infeasible!\n" << std::endl;
        throw qol::Exception("Infeasible solution!");
      }

    }

      std::string outtimefile = "times/" + probModel->getName()+"_times.csv";
      std::ofstream outtime;
		  outtime.open(outtimefile, std::ios_base::app);
		  outtime << sh.WINDOW_SIZE << "," << timer.elapsedWallTime()
              << "," << orig_obj << "," << obj << std::endl;
      std::cout << "pass time written to " << outtimefile << ". ending run!"
                << std::endl << std::endl;

  }

    sh.WINDOW_SIZE++;
    dir_indicator=true;
  }

  return obj;
}

double SinglePSolver::slideWindow(std::vector<std::vector<std::vector<double> > > &Y_sol){

  qol::CpuTimer timer;

  const int t_max = probModel->getNPeriod();
  const double eps = 1e-6;

  bool early_finish = false;

  std::cout << "\n*** running slideWindow() ***\n"<< std::endl;

  std::cout << "Window size: " << sh.WINDOW_SIZE << std::endl;
  std::vector<double> objectives;
  std::vector<double> times;

  double obj = 0.0;
  double old_obj = -1e9;

	int converge_count = 0;

  for (size_t i=0;i<sh.NUM_PASSES && !early_finish;++i){

    std::cout << "\n*** running backward pass "<< i+1 << " ***" << std::endl;

    for (size_t t=t_max-sh.WINDOW_SIZE-1; t>0 && !early_finish; --t){
      obj = improveWindow(Y_sol, t, i, sh.WINDOW_SIZE, false);
      std::cout << "\nImproved objective: " << obj << std::endl;

      if (sh.OBJ_UB > 0)
        std::cout << "\nObjective upper bound: " << sh.OBJ_UB << std::endl;

      objectives.push_back(obj);
      times.push_back(timer.elapsedWallTime());

      if (sh.RECORD_DATA){
        std::ofstream outfile;
        outfile.open ("timeobj.csv");
        for (size_t i=0; i<objectives.size(); ++i)
          outfile << times[i] << "," << objectives[i] << std::endl;
        outfile.close();
      }

      std::cout << std::endl << "WINDOW SEARCH TOTAL TIME: " << timer.elapsedSeconds()
                << " sec CPU / " << timer.elapsedWallTime() << " sec wall "
                << std::endl;

      if (sh.WALL_SEARCH_TIME > 0)
        std::cout << "\nTime limit (wall): " << sh.WALL_SEARCH_TIME << std::endl;

      if (sh.CPU_SEARCH_TIME > 0)
        std::cout << "\nTime limit (CPU): " << sh.CPU_SEARCH_TIME << std::endl;

      Sol_Int test_sol;
      packSol(Y_sol, test_sol);

      test_sol.nT = t_max;

      std::cout << "\nChecking solution from solver...\n" << std::endl;

      bool test_error = verify((*probModel), test_sol);

      std::cout << "Solution found by solver was ";
      if (!test_error)
        std::cout << "feasible!\n" << std::endl;
      else{
        std::cout << "infeasible!\n" << std::endl;
        throw qol::Exception("Infeasible solution!");
      }

      if (obj > sh.OBJ_UB + eps && sh.OBJ_UB > 0){
        std::cout << "\nobjective: " << obj << " over upper bound of "
                  << sh.OBJ_UB << ", early finish!!" << std::endl;
        early_finish=true;
      }

      if (sh.WALL_SEARCH_TIME > 0 && timer.elapsedWallTime() > sh.WALL_SEARCH_TIME){
        std::cout << std::endl << "wall time limit of " << sh.WALL_SEARCH_TIME
                  << " exceeded, early finish!!" << std::endl;
        early_finish=true;
      }
      if (sh.CPU_SEARCH_TIME > 0 && timer.elapsedSeconds() > sh.CPU_SEARCH_TIME){
        std::cout << std::endl << "CPU time limit of " << sh.CPU_SEARCH_TIME
                  << " exceeded, early finish!!" << std::endl;
        early_finish=true;
      }
    }

    if (obj - old_obj > eps){
      old_obj = obj;
			converge_count = 0;
    }
    else{
      converge_count++;
      std::cout << "********** NO IMPROVEMENT IN " << converge_count
                << " FULL PASSES **********\n"<< std::endl;

		}

		if (sh.CONVERGE_STOP && converge_count >= sh.CONVERGE_STOP) {
      std::cout << "********** CONVERGE LIMIT REACHED "
                << " ENDING RUN **********\n"<< std::endl;
      early_finish = true;
    }

    if (sh.RECORD_PASS_TIME){
      std::string outtimefile = "times/" + probModel->getName()+"_times.csv";
      std::ofstream outtime;
		  outtime.open(outtimefile, std::ios_base::app);
		  outtime << sh.WINDOW_SIZE << "," << timer.elapsedWallTime() << std::endl;
      std::cout << "pass time written to " << outtimefile << ". ending run!"
                << std::endl << std::endl;
      early_finish = true;
    }

    if (!early_finish)
      std::cout << "*** running foward pass " << i+1 << " ***" << std::endl;

    for (size_t t=0; t<t_max-sh.WINDOW_SIZE && !early_finish; ++t){
      obj = improveWindow(Y_sol, t, i, sh.WINDOW_SIZE, true);
      std::cout << "\nImproved objective: " << obj << std::endl;

      if (sh.OBJ_UB > 0)
        std::cout << "\nObjective upper bound: " << sh.OBJ_UB << std::endl;

      objectives.push_back(obj);
      times.push_back(timer.elapsedWallTime());

      if (sh.RECORD_DATA){
        std::ofstream outfile;
        outfile.open ("timeobj.csv");
        for (size_t i=0; i<objectives.size(); ++i)
          outfile << times[i] << "," << objectives[i] << std::endl;
        outfile.close();
      }

      std::cout << std::endl << "WINDOW SEARCH TOTAL TIME: " << timer.elapsedSeconds()
                << " sec CPU / " << timer.elapsedWallTime() << " sec wall "
                << std::endl;

      if (sh.WALL_SEARCH_TIME > 0)
        std::cout << "\nTime limit (wall): " << sh.WALL_SEARCH_TIME << std::endl;

      if (sh.CPU_SEARCH_TIME > 0)
        std::cout << "\nTime limit (CPU): " << sh.CPU_SEARCH_TIME << std::endl;


      std::cout << "\nChecking solution from solver...\n" << std::endl;

      Sol_Int test_sol;
      packSol(Y_sol, test_sol);

      test_sol.nT = t_max;

      bool test_error = verify((*probModel), test_sol);

      std::cout << "Solution found by solver was ";
      if (!test_error)
        std::cout << "feasible!\n" << std::endl;
      else{
        std::cout << "infeasible!\n" << std::endl;
        throw qol::Exception("Infeasible solution!");
      }

      if (obj > sh.OBJ_UB + eps && sh.OBJ_UB > 0){
        std::cout << "\nobjective: " << obj << " over upper bound of "
                  << sh.OBJ_UB << ", early finish!!" << std::endl;
        early_finish=true;
      }

      if (sh.WALL_SEARCH_TIME > 0 && timer.elapsedWallTime() > sh.WALL_SEARCH_TIME){
        std::cout << std::endl << "wall time limit of " << sh.WALL_SEARCH_TIME
                  << " exceeded, early finish!!" << std::endl;
        early_finish=true;
      }

      if (sh.CPU_SEARCH_TIME > 0 && timer.elapsedSeconds() > sh.CPU_SEARCH_TIME){
        std::cout << std::endl << "CPU time limit of " << sh.CPU_SEARCH_TIME
                  << " exceeded, early finish!!" << std::endl;
        early_finish=true;
      }
    }

    if (!early_finish){
      if (obj - old_obj > eps){
        old_obj = obj;
  			converge_count = 0;
      }
      else{
        converge_count++;
        std::cout << "********** NO IMPROVEMENT IN " << converge_count
                  << " FULL PASSES **********\n"<< std::endl;

  		}

  		if (sh.CONVERGE_STOP && converge_count >= sh.CONVERGE_STOP) {
        std::cout << "********** CONVERGE LIMIT REACHED "
                  << " ENDING RUN **********\n"<< std::endl;
        early_finish = true;
      }
    }

  }

  if (sh.RECORD_DATA)
    std::cout << "\nslideWindow() complete! Data available in timeobj.csv" << std::endl;
  else{
    std::cout << "\nslideWindow() complete. time, objective values:" << std::endl;
    for (size_t i=0; i<objectives.size(); ++i)
      std::cout << times[i] << ", " << objectives[i] << std::endl;
  }
  std::cout << std::endl;

  return obj;
}

double SinglePSolver::improveWindow(std::vector<std::vector<std::vector<double> > > &Y_sol,
                                    int time, int pass_num, int window_size,
                                    bool forward_pass){
  // set up some parameters
  const int nB = probModel->getNBlock();
  const int t_max = probModel->getNPeriod();
  const int d_max = probModel->getnDestination();
  const int r_max = probModel->getnResources();
  const double rate = probModel->getDiscountRate();
  std::vector<Block> * blocks=probModel->getBlock();

  bool last_window = false;
  double t_obj = 0.0;

  std::vector<std::vector<double> > Xt(nB,
      std::vector<double>(window_size,0.0)); // X vector over time horizon

  std::vector<bool> in_window(nB, false); // freed blocks within time points
  std::vector<bool> freed(nB, false);
  std::vector<bool> no_process(nB, false); // blocks with value below processing threshold

  std::cout << "Fixing blocks ";

  if (time > 0)
    std::cout << "before " << time << " and ";

  std::cout << "after " << time+window_size-1;

  if (forward_pass)
    std::cout << " - forward pass " << pass_num;
  else
    std::cout << " - backward pass " << pass_num;

  std::cout << std::endl;

  std::cout << "getting blocks in window.. " << std::flush;

  int b_count = getBlocksInWindow(Y_sol, Xt, in_window, no_process, time);

  std::cout << "done!" << std::endl;

  int count_verify = 0;

  for (size_t b=0;b<nB;++b){
    if (in_window[b])
      count_verify++;
  }

  std::cout << b_count << " blocks freed (verified " << count_verify << ")" << std::endl;

  int old_b_count = b_count;

  if (sh.REDUCED_BLOCKS > 0){
    b_count = getReducedBlocks_NoBoundary(Y_sol, in_window, freed, time, b_count);
  }else{
    freed=in_window;
  }

  count_verify = 0;

  for (size_t b=0;b<nB;++b){
    if (freed[b])
      count_verify++;
  }

  std::cout << b_count << " blocks freed, reduced from "
  << old_b_count << " (verified " << count_verify << ")" << std::endl;

//  int b_count = 0; // freed block count
//  for (size_t b=0; b<nB; ++b){
//    for (size_t t=time; t<time+sh.WINDOW_SIZE;++t){ // iterate only within window
//      for (size_t d=0; d<d_max; ++d){
//        if (Y_sol[b][t][d] > 0.0){
//          //if (sh.ONLY_NEGATIVE){
//          //  double profit=(*blocks)[b].getProfit(d);
//          //  if (profit > 0){
//                Xt[b][t-time] = 1.0;
//          //    }
//          //}
//          freed[b] = true;
//          b_count++;
//        }
//      }
//    }
//  }

  //for (int b=0;b<Xt.size();++b){
  //    std::cout<<"block: " << b << " -";
  //  for (int i=0;i<Xt[b].size();++i){
  //    std::cout << " " << Xt[b][i];
  //  }
  //  std::cout << std::endl;
  //}


  t_obj = MIP_window(in_window, freed, time, time+window_size, Y_sol, Xt, last_window, b_count);

  std::cout << "repairing objective.." << std::endl;

  Sol_Int temp_sol;

  packSol(Y_sol, temp_sol);

  t_obj = getTrueObj(temp_sol);

  return t_obj;

}

// TODO: TRY NOT USING BOUNDARY, JUST CALCULATE CONE VALUES
// AND ADD HIGH VALUE CONES
// time is absolute time
int SinglePSolver::getReducedBlocks(
                   const std::vector<std::vector<std::vector<double> > >Y_sol,
                   std::vector<bool> freed, std::vector<bool> &reduced_freed, int time, int blocks_in_window){

  const int nB = probModel->getNBlock();
  const int d_max = probModel->getnDestination();
  const int t_max = probModel->getNPeriod();
  const double rate = probModel->getDiscountRate();

  int b_count = 0;
  std::vector<bool> on_L_boundary (nB, false); // on t,t-1
  std::vector<bool> on_R_boundary (nB, false); // on t+WINDOW_SIZE,t+WINDOW_SIZE+1
  std::vector<std::vector<bool> > inCone (nB, std::vector<bool> (nB, false));
  std::vector<double> cone_value_L (nB, 0.0);
  std::vector<double> cone_value_R (nB, 0.0);

  int extremeL = 0;
  int extremeR = 0;

  int L_B_size = 0;
  int R_B_size = 0;

  std::cout << "finding boundaries.. " << std::flush;
#   pragma omp parallel for
  for (size_t b=0; b<nB; ++b){
    if (freed[b]){ // if block is in current window
      const Node* curr = probModel->graph.getNode(b);
      bool found = false;

      // // check for extreme boundaries
      // if (curr->getInDegree() == 0){
      //   on_L_boundary[b] = true;
      //   L_B_size++;
      //   found = true;
      //   extremeL++;
      // }
      // if (curr->getOutDegree() == 0){
      //   on_R_boundary[b] = true;
      //   R_B_size++;
      //   found = true;
      //   extremeR++;
      // }

      // check for left boundary
      // block is on boundary if it has predecessors in the previous time window

      //TODO FIND BOUNDARIES PROPERLY!

      if (time > 0 && curr->getInDegree() > 0){
        for (size_t i=0; i<curr->getInDegree() && !found; ++i){
          int prev = curr->getInArc(i)->getSrcID();
          for (size_t d=0; d<d_max && !found; ++d){
            if (Y_sol[prev][time-1][d] > 0){
              on_L_boundary[b] = true;
              L_B_size++;
              found = true;
            }
          }
        }
      }
      else {
        if (curr->getInDegree() == 0){
          on_L_boundary[b] = true;
          L_B_size++;
          found = true;
          if (time == 0)
            extremeL++;
        }
      }
      // check for right boundary
      // block is on boundary if one or more of its successors are not in the current
      // time window

      found = false;

      if (time < t_max && curr->getOutDegree() > 0){
        for (size_t i=0; i<curr->getOutDegree() && !found; ++i){
          int next = curr->getOutArc(i)->getSrcID();
          for (size_t d=0; d<d_max && !found; ++d){
            if (Y_sol[next][time+1][d] > 0){
              on_R_boundary[b] = true;
              R_B_size++;
              found = true;
            }
          }
        }
      }
      else {
        if (curr->getOutDegree() == 0){
          on_R_boundary[b] = true;
          R_B_size++;
          found = true;
          if (time == t_max)
            extremeR++;
        }
      }
    }
  }

  std::cout << "done!" << std::endl;

  std::cout << "validating R and L boundary blocks.. ";

  int duplicate_count=0;

  for (size_t b=0; b<nB; ++b){
    if (on_L_boundary[b] > 0 && on_R_boundary[b] > 0)
      duplicate_count++;
      //std::cout << "ERROR! block " << b << " on both boundaries!" << std::endl;
  }

  std::cout << "done!" << std::endl << duplicate_count
            << " blocks on both boundaries" << std::endl;

  std::cout
    << "number of blocks on L boundary: " << L_B_size << ", extremes: " << extremeL << std::endl
    << "number of blocks on R boundary: " << R_B_size << ", extremes: " << extremeR << std::endl << std::endl
    << "getting left hand cones.. " << std::flush;

  // get cone values for left hand boundary blocks within the same half window
#   pragma omp parallel for
  for (size_t b=0; b<nB; ++b){
    if (on_L_boundary[b]){
      std::vector<int> stack;
      stack.push_back(b);
      inCone[b][b] = true;
      while (!stack.empty()){
        int blockID = stack.back();
        const Node* curr = probModel->graph.getNode(blockID);
        const Block & block=probModel->getBlock(blockID);
        stack.pop_back();
        double max_profit = block.getProfit(0)/pow(1+rate,time);
        for (int d=1;d<d_max;++d){
          max_profit = std::max(max_profit, block.getProfit(d)/pow(1+rate,time));
        }
        // subtract profit to allow for ordering using priority queue
        cone_value_L[b] += max_profit;

        // iterate over successors
        for (size_t i=0; i < curr->getOutDegree(); ++i){
          int next = curr->getOutArc(i)->getTgtID();
          if (next < nB && next >= 0){
            bool same_period = false;
            // iterate over half the window size (ceiling)
            for (size_t t=time;t<time+ceil(double(sh.WINDOW_SIZE)/2) && !same_period; ++t){
              for (size_t d=0;d<d_max && !same_period;++d){
                if (Y_sol[next][t][d] > 0.0){
                  same_period = true;
                }
              }
            }
            if (!inCone[b][next] && same_period){
              inCone[b][next] = true;
              stack.push_back(next);
            }
          }
        }
      }
    }
  }

  std::cout << "done!" << std::endl << "getting right hand cones.. " << std::flush;

  // get cone values for right hand boundary blocks within same period
#   pragma omp parallel for
  for (size_t b=0; b<nB; ++b){
    if (on_R_boundary[b]){
      std::vector<int> stack;
      stack.push_back(b);
      inCone[b][b] = true;
      while (!stack.empty()){
        int blockID = stack.back();
        const Node* curr = probModel->graph.getNode(blockID);
        const Block & block=probModel->getBlock(blockID);
        stack.pop_back();
        double max_profit = block.getProfit(0)/pow(1+rate,time+sh.WINDOW_SIZE-1);
        for (int d=1;d<d_max;++d){
          max_profit = std::max(max_profit,
              block.getProfit(d)/pow(1+rate,time+sh.WINDOW_SIZE-1));
        }
        cone_value_R[b] += max_profit;

        // iterate over successors
        for (size_t i=0; i < curr->getInDegree(); ++i){
          int prev = curr->getInArc(i)->getSrcID();
          bool same_period = false;
          for (size_t t=time+(sh.WINDOW_SIZE/2);t<time+sh.WINDOW_SIZE;++t){
            for (size_t d=0;d<d_max && !same_period;++d){
              if (Y_sol[prev][t][d] > 0.0){
                same_period = true;
              }
            }
          }
          if (!inCone[b][prev] && same_period){
            inCone[b][prev] = true;
            stack.push_back(prev);
          }
        }
      }
    }
  }

  std::cout << "done!" << std::endl << "ordering cones.. " << std::flush;

  std::priority_queue<std::pair<double, int> > q_L;
  std::priority_queue<std::pair<double, int> > q_R;

  double L_list_min = 9e10;
  double R_list_min = 9e10;

  for (size_t b = 0; b<nB; ++b){
    if (on_L_boundary[b]){
      if (cone_value_L[b] > L_list_min || b < sh.REDUCED_BLOCKS*nB){
        q_L.push(std::pair<double, int> (cone_value_L[b],b));
        if (cone_value_L[b] < L_list_min)
          L_list_min = cone_value_L[b];
      }
    }
    if (on_R_boundary[b]){
      if (cone_value_R[b] > R_list_min || b < sh.REDUCED_BLOCKS*nB){
        q_R.push(std::pair<double, int> (cone_value_R[b],b));
        if (cone_value_R[b] < R_list_min)
          R_list_min = cone_value_R[b];
      }
    }
  }
  std::cout << "done!" << std::endl;

  // take cones from L and R boundary in an alternating fashion until
  // total number of blocks freed are sh.REDUCED_BLOCKS*blocks_in_window
  // where blocks_in_window represents the original number of blocks freed

  int count_verify = 0;
  for (size_t b=0;b<nB;++b){
    if (reduced_freed[b])
      count_verify++;
  }

  std::cout << "done!" << std::endl << "count_verify: " << count_verify << std::endl;

  bool use_L = true;
  int L_cones_used = 0;
  int R_cones_used = 0;

  while (b_count < blocks_in_window*sh.REDUCED_BLOCKS && (q_L.size() > 0 || q_R.size() > 0)){
    int curr;
    if (use_L && q_L.size() > 0){
      curr = q_L.top().second;
      L_cones_used++;
      q_L.pop();
    }else{
      if (q_R.size() > 0){
        curr = q_R.top().second;
        R_cones_used++;
        q_R.pop();
      }
    }

    for (size_t b=0;b<nB;++b){
      if (!reduced_freed[b] && inCone[curr][b]){
        reduced_freed[b] = true;
        b_count++;
      }
    }

    use_L = !use_L;
  }

  std::cout << L_cones_used << " L cones added" << std::endl
            << R_cones_used << " R cones added" << std::endl;

  // std::cout << q_L.size()
  //           << " cones in L boundary set, printing top 10 L cones: " << std::endl;
  //
  // for (int i=0;i<10 && q_L.size() > 0;++i){
  //   int cone_size = 0;
  //   for (int b=0;b<inCone[q_L.top().second].size();++b){
  //     if(inCone[q_L.top().second][b])
  //       cone_size++;
  //   }
  //   std::cout << "block: " << q_L.top().second << ", value: "
  //             << q_L.top().first << ", cone_size: " << cone_size << std::endl;
  //  q_L.pop();
  // }
  //
  // std::cout << q_R.size()
  //           << " cones in R boundary set, printing top 10 R cones: " << std::endl;
  //
  // for (int i=0;i<10 && q_R.size() > 0;++i){
  //   int cone_size = 0;
  //   for (int b=0;b<inCone[q_R.top().second].size();++b){
  //     if(inCone[q_R.top().second][b])
  //       cone_size++;
  //   }
  //   std::cout << "block: " << q_R.top().second << ", value: "
  //             << q_R.top().first << ", cone_size: " << cone_size << std::endl;
  //  q_R.pop();
  // }

  count_verify = 0;
  for (size_t b=0;b<nB;++b){
    if (reduced_freed[b])
      count_verify++;
  }

  std::cout << "count_verify: " << count_verify << std::endl;


  return b_count;
}

int SinglePSolver::getReducedBlocks_NoBoundary(
                   const std::vector<std::vector<std::vector<double> > >Y_sol,
                   const std::vector<bool> &in_window, std::vector<bool> &freed, int time, int blocks_in_window){

  const int nB = probModel->getNBlock();
  const int d_max = probModel->getnDestination();
  const int t_max = probModel->getNPeriod();
  const double rate = probModel->getDiscountRate();

  int b_count = 0;
  std::vector<bool> on_L_boundary (nB, false); // on t
  std::vector<bool> on_R_boundary (nB, false); // on t+WINDOW_SIZE-1
  std::vector<std::vector<bool> > inCone (nB, std::vector<bool> (nB, false));
  std::vector<double> cone_value_L (nB, 0.0);
  std::vector<double> cone_value_R (nB, 0.0);

  int L_B_size = 0;
  int R_B_size = 0;

  //TODO: for some reason the boundary sets arent adding up!
  //      when finding all boundary blocks, mark off which ones have been found
  //      and compare against master freed list, see which one is missing!

  std::cout << "finding boundaries.. " << std::flush;
#   pragma omp parallel for
  for (size_t b=0; b<nB; ++b){
    if (in_window[b]){ // if block is in current window
      bool found = false;

      for (size_t d=0; d<d_max && !found; ++d){
        if(Y_sol[b][time][d] > 0.0){
          on_L_boundary[b] = true;
          L_B_size++;
          found = true;
        }
        if(Y_sol[b][time+sh.WINDOW_SIZE-1][d] > 0.0 && !found){
          on_R_boundary[b] = true;
          R_B_size++;
          found = true;
        }
      }
    }
  }

  std::cout << "done!" << std::endl;

  std::cout << "validating R and L boundary blocks.. ";

  int duplicate_count=0;

  int count_verify_L = 0;
  int count_verify_R = 0;

  for (size_t b=0; b<nB; ++b){
    if (on_L_boundary[b] > 0 && on_R_boundary[b] > 0)
      duplicate_count++;
      //std::cout << "ERROR! block " << b << " on both boundaries!" << std::endl;
    if (on_L_boundary[b])
      count_verify_L++;
    if (on_R_boundary[b])
      count_verify_R++;
  }

  std::cout << "done!" << std::endl << duplicate_count
            << " blocks on both boundaries" << std::endl;

  std::cout
    << "number of blocks on L boundary: " << L_B_size << " (" << count_verify_L
    << " verified)" << std::endl
    << "number of blocks on R boundary: " << R_B_size << " (" << count_verify_R
    << " verified)"<< std::endl
    << "total sum of L and R boundaries: " << (L_B_size + R_B_size) << " ("
    << (count_verify_R + count_verify_L)   << " verified)" << std::endl;

  std::vector<bool> not_on_boundary (nB, false);
  int nobound_count = 0;

  for (size_t b=0;b<nB;++b){
    bool found_once = false;
    if (in_window[b]){
      if (on_L_boundary[b])
        found_once = true;
      if (on_R_boundary[b]){
        if (found_once){
          std::cout << "block " << b << " on both boundaries!" << std::endl;
        } else {
          found_once = true;
        }
      }
      if (!found_once){
        not_on_boundary[b]=true;
        nobound_count++;
      }
    }
  }

  std::cout << nobound_count << " blocks not on either boundary" << std::endl;

  std::cout
    << "getting left hand cones.. " << std::flush;

  // get cone values for left hand boundary blocks within the same half window
  // negative profit because want low valued cones from right
#   pragma omp parallel for
  for (size_t b=0; b<nB; ++b){
    if (on_L_boundary[b]){
      std::vector<int> stack;
      stack.push_back(b);
      inCone[b][b] = true;
      while (!stack.empty()){
        int blockID = stack.back();
        const Node* curr = probModel->graph.getNode(blockID);
        const Block & block=probModel->getBlock(blockID);
        stack.pop_back();
        double max_profit = block.getProfit(0)/pow(1+rate,time);
        for (int d=1;d<d_max;++d){
          max_profit = std::max(max_profit, block.getProfit(d)/pow(1+rate,time));
        }
        // subtract profit to allow for ordering using priority queue
        cone_value_L[b] -= max_profit;

        // iterate over successors
        for (size_t i=0; i < curr->getOutDegree(); ++i){
          int next = curr->getOutArc(i)->getTgtID();
          if (next < nB && next >= 0){
            bool same_period = false;
            // iterate over half the window size (ceiling)
            for (size_t t=time;t<time+ceil(double(sh.WINDOW_SIZE)/2) && !same_period; ++t){
              for (size_t d=0;d<d_max && !same_period;++d){
                if (Y_sol[next][t][d] > 0.0){
                  same_period = true;
                }
              }
            }
            if (!inCone[b][next] && same_period){
              inCone[b][next] = true;
              stack.push_back(next);
            }
          }
        }
      }
    }
  }

  std::cout << "done!" << std::endl << "getting right hand cones.. " << std::flush;

  // get cone values for right hand boundary blocks within same period
#   pragma omp parallel for
  for (size_t b=0; b<nB; ++b){
    if (on_R_boundary[b]){
      std::vector<int> stack;
      stack.push_back(b);
      inCone[b][b] = true;
      while (!stack.empty()){
        int blockID = stack.back();
        const Node* curr = probModel->graph.getNode(blockID);
        const Block & block=probModel->getBlock(blockID);
        stack.pop_back();
        double max_profit = block.getProfit(0)/pow(1+rate,time+sh.WINDOW_SIZE-1);
        for (int d=1;d<d_max;++d){
          max_profit = std::max(max_profit,
              block.getProfit(d)/pow(1+rate,time+sh.WINDOW_SIZE-1));
        }
        cone_value_R[b] += max_profit;

        // iterate over successors
        for (size_t i=0; i < curr->getInDegree(); ++i){
          int prev = curr->getInArc(i)->getSrcID();
          bool same_period = false;
          for (size_t t=time+(sh.WINDOW_SIZE/2);t<time+sh.WINDOW_SIZE;++t){
            for (size_t d=0;d<d_max && !same_period;++d){
              if (Y_sol[prev][t][d] > 0.0){
                same_period = true;
              }
            }
          }
          if (!inCone[b][prev] && same_period){
            inCone[b][prev] = true;
            stack.push_back(prev);
          }
        }
      }
    }
  }

  std::cout << "done!" << std::endl << "ordering cones.. " << std::flush;

  std::priority_queue<std::pair<double, int> > q_L;
  std::priority_queue<std::pair<double, int> > q_R;

  double L_list_min = 9e10;
  double R_list_min = 9e10;

  for (size_t b = 0; b<nB; ++b){
    if (on_L_boundary[b]){
      if (cone_value_L[b] > L_list_min || q_L.size() < ceil(sh.REDUCED_BLOCKS*blocks_in_window)){
        q_L.push(std::pair<double, int> (cone_value_L[b],b));
        if (cone_value_L[b] < L_list_min)
          L_list_min = cone_value_L[b];
      }
    }
    if (on_R_boundary[b]){
      if (cone_value_R[b] > R_list_min || q_R.size() < ceil(sh.REDUCED_BLOCKS*blocks_in_window)){
        q_R.push(std::pair<double, int> (cone_value_R[b],b));
        if (cone_value_R[b] < R_list_min)
          R_list_min = cone_value_R[b];
      }
    }
  }
  std::cout << "done!" << std::endl;

  // take cones from L and R boundary in an alternating fashion until
  // total number of blocks freed are sh.REDUCED_BLOCKS*blocks_in_window
  // where blocks_in_window represents the original number of blocks freed

  int count_verify = 0;
  for (size_t b=0;b<nB;++b){
    if (freed[b])
      count_verify++;
  }

  std::cout << "done!" << std::endl << "count_verify: " << count_verify << std::endl;

  bool use_L = true;
  int L_cones_used = 0;
  int R_cones_used = 0;

  std::cout << blocks_in_window << " blocks in window" << std::endl
            << "adding " << ceil(blocks_in_window*sh.REDUCED_BLOCKS)
            << " blocks" << std::endl << "q_L.size() = " << q_L.size() << std::endl
            << "q_R.size() = " << q_R.size() << std::endl;

  while (b_count < ceil(blocks_in_window*sh.REDUCED_BLOCKS) && (q_L.size() > 0 || q_R.size() > 0)){
    int curr;
    if (use_L && q_L.size() > 0){
      curr = q_L.top().second;
      L_cones_used++;
      q_L.pop();
    }
    else{
      if (q_R.size() > 0){
        curr = q_R.top().second;
        R_cones_used++;
        q_R.pop();
      }
    }

    for (size_t b=0;b<nB;++b){
      if (!freed[b] && inCone[curr][b]){
        freed[b] = true;
        b_count++;
      }
    }

    use_L = !use_L;
  }

  std::cout << L_cones_used << " L cones added" << std::endl
            << R_cones_used << " R cones added" << std::endl;

  // std::cout << q_L.size()
  //           << " cones in L boundary set, printing top 10 L cones: " << std::endl;
  //
  // for (int i=0;i<10 && q_L.size() > 0;++i){
  //   int cone_size = 0;
  //   for (int b=0;b<inCone[q_L.top().second].size();++b){
  //     if(inCone[q_L.top().second][b])
  //       cone_size++;
  //   }
  //   std::cout << "block: " << q_L.top().second << ", value: "
  //             << q_L.top().first << ", cone_size: " << cone_size << std::endl;
  //  q_L.pop();
  // }
  //
  // std::cout << q_R.size()
  //           << " cones in R boundary set, printing top 10 R cones: " << std::endl;
  //
  // for (int i=0;i<10 && q_R.size() > 0;++i){
  //   int cone_size = 0;
  //   for (int b=0;b<inCone[q_R.top().second].size();++b){
  //     if(inCone[q_R.top().second][b])
  //       cone_size++;
  //   }
  //   std::cout << "block: " << q_R.top().second << ", value: "
  //             << q_R.top().first << ", cone_size: " << cone_size << std::endl;
  //  q_R.pop();
  // }

  count_verify = 0;
  for (size_t b=0;b<nB;++b){
    if (freed[b])
      count_verify++;
  }

  std::cout << "count_verify: " << count_verify << std::endl;


  return b_count;
}

// time is absolute time
int SinglePSolver::getBlocksInWindow(
                  const std::vector<std::vector<std::vector<double> > > &Y_sol,
                  std::vector<std::vector<double> > &Xt,
                  std::vector<bool> &in_window, std::vector<bool> &no_process, int time){

  const int nB = probModel->getNBlock();
  const int d_max = probModel->getnDestination();



  int b_count = 0; // freed block count
  for (size_t b=0; b<nB; ++b){
    for (size_t t=time; t<time+sh.WINDOW_SIZE;++t){ // iterate only within window
      bool found = false;
      for (size_t d=0; d<d_max && !found; ++d){
        if (Y_sol[b][t][d] > 0.0){
          //if (sh.ONLY_NEGATIVE){
          //  double profit=(*blocks)[b].getProfit(d);
          //  if (profit > 0){
                Xt[b][t-time] = 1.0;
          //    }
          //}
          in_window[b] = true;
          found = true;
          b_count++;
        }
      }
    }
  }

  return b_count;
}

// TODO: work out what is going on with the constraints and stuff,
//       its something to do with the fact that you can have a difference
//       between the freed set and the in_window set, but not sure what
//       it is!

double SinglePSolver::MIP_window(const std::vector<bool> &in_window,
                      const std::vector<bool> &freed,
                      int mod_tmin, int mod_tmax,
                      std::vector<std::vector<std::vector<double> > > &Y_sol,
                      const std::vector<std::vector<double> > &X,
                      bool last_window, int b_count){

  qol::CpuTimer timer;

  // set up some parameters
  const int nB = probModel->getNBlock();
  const int t_max = probModel->getNPeriod();
  const int d_max = probModel->getnDestination();
  const int r_max = probModel->getnResources();
  const double rate = probModel->getDiscountRate();
  std::vector<Block> * blocks=probModel->getBlock();

  std::string lpFile="";// ,solnFile = probModel->getName()+".sol";

  qol::MIPSolver *mipPtrW = 0;

  qol::Parameters param;
  param.setParamVal(qol::VERBOSITY,1);
  if (sh.WINDOW_SEARCH_TIME > 0)
    param.setParamVal(qol::TIMELIMIT,sh.WINDOW_SEARCH_TIME);
  //param.setParamVal(qol::RELGAP, 0.02);

  mipPtrW = new qol::CplexFormulation();

  qol::MIPSolver &model = *mipPtrW;

  model.setParameters(param);

  std::vector<std::vector<qol::Variable> > Xvar(nB);
  std::vector<std::vector<std::vector<qol::Variable> > > Yvar(nB);

  std::cout << "Setting up extended model ... " << std::endl;

  //int time = mod_tmax - mod_tmin;
  double obj = 0.0;

  // create variables and set them to binary
  // only allow positive values for blocks that exist
  for (size_t b=0; b<nB; ++b){
    Xvar[b].resize(sh.WINDOW_SIZE);
    Yvar[b].resize(sh.WINDOW_SIZE);

    // TODO: get profit and adjust it for current time period maybe
    //       by switching the order between d and t iterators otherwise just
    //       compute profit on the fly - DONE (computed on fly)

    for (size_t t=0; t<sh.WINDOW_SIZE; t++){
      std::ostringstream Xvar_label;
      Xvar_label << "x_" << b << "_" << t << std::endl;

      if (in_window[b]){
        if (freed[b]){
          //if(X[b][t] == 1.0)
          //  Xvar[b][t] = model.addVar(1,1,0,qol::Variable::BINARY, Xvar_label.str());
          //else
          Xvar[b][t] = model.addVar(0,1,0,qol::Variable::BINARY, Xvar_label.str());
        }
        else{
          Xvar[b][t] = model.addVar(X[b][mod_tmin+t],X[b][mod_tmin+t],0,qol::Variable::BINARY, Xvar_label.str());
        }
      }
      else
        Xvar[b][t] = model.addVar(0,0,0,qol::Variable::BINARY, Xvar_label.str());

      //std::cout << "adding xvar.. " << std::flush;
      //if (WARM_START)
      //  model.setPrimalStart(Xvar[b][t], X[b][t]);
      //std::cout << "done!" << std::endl;


      Yvar[b][t].resize(d_max);
      for (size_t d=0;d<d_max;++d){
        std::ostringstream Yvar_label;
        Yvar_label << "y_" << b << "_" << t << "_" << d << std::endl;

        if (in_window[b]){
          double block_profit = (*blocks)[b].getProfit(d);
          double profit = block_profit/pow(1+rate,t);
          if (freed[b]){
            Yvar[b][t][d] = model.addVar(0,1,-profit,qol::Variable::CONTINUOUS, Yvar_label.str());
          }
          else{
            Yvar[b][t][d] = model.addVar(Y_sol[b][mod_tmin+t][d],Y_sol[b][mod_tmin+t][d],-profit,qol::Variable::CONTINUOUS, Yvar_label.str());
          }
        }
        else
          Yvar[b][t][d] = model.addVar(0,0,0,qol::Variable::CONTINUOUS, Yvar_label.str());

        //std::cout << "adding yvar.. " << std::flush;
        //if (WARM_START)
        //  model.setPrimalStart(Yvar[b][t][d], Y_sol[b][t+mod_tmin][d]);
        //std::cout << "done!" << std::endl;
      }
    }
  }

  if (sh.WARM_START){
    for (size_t b=0;b<nB;++b){
      for (size_t t=0;t<sh.WINDOW_SIZE;++t){
        model.setPrimalStart(Xvar[b][t],X[b][mod_tmin+t]);
        for (size_t d=0;d<d_max;++d){
          model.setPrimalStart(Yvar[b][t][d],Y_sol[b][mod_tmin+t][d]);
        }
      }
    }
  }


  // constraints
  // precedence (7)
  for (size_t a=0; a<nB; ++a){
    if (!in_window[a])
      continue; // block a not used
    std::vector<int> * pred = (*blocks)[a].getPreds();
    int n = (*blocks)[a].getNumPred();
    for (size_t p=0; p<n; ++p){
      int b = (*pred)[p];
      if (!in_window[b])
        continue;
      for (size_t t=0; t<sh.WINDOW_SIZE; t++){
        model.addConstraint(Xvar[a][t] <= Xvar[b][t]).setName("x_%d_%d_leq_x_%d_%d",a,t,b,t);
      }
    }
  }

  // SumDest (8)
  for (size_t b=0; b<nB; ++b){
    if(!freed[b])
      continue;
    for (size_t t=0; t<sh.WINDOW_SIZE; ++t){
      qol::Expression sumY;
      for (size_t d=0; d<d_max; ++d){
        sumY += Yvar[b][t][d];
      }
      qol::Expression rhs = Xvar[b][t];
      if (t > 0){
        rhs -= Xvar[b][t-1];
      }
      model.addConstraint( sumY == rhs ).setName("dest_%d_%d",b,t);
    }
  }

  // block remains done (9)
  for (size_t b = 0; b<nB; ++b){
    for (size_t t = 0; t<sh.WINDOW_SIZE-1; ++t){
      std::ostringstream remain_label;
      remain_label << "xt_" << b << "_" << t << std::endl;
      model.addConstraint( Xvar[b][t] <= Xvar[b][t+1] ).setName("xt_%d_%d",b,t);
    }
  }

  // last window constraint to ensure every block completes
  for(size_t b=0; b<nB && !last_window; b++){
    if(freed[b]){
      qol::Expression xbt_minus1 = Xvar[b][sh.WINDOW_SIZE-1];
      model.addConstraint( xbt_minus1 == 1).setName("lxw_%d",b);
    }
  }

  // Resource constraints (10)
  for(size_t r=0; r<r_max; r++){
    for (size_t t=0; t<sh.WINDOW_SIZE; ++t){
      char cType = probModel->getResConstrType(r, t);
      qol::Expression expr;
      for(size_t b=0; b<nB; b++){
        if (!freed[b])
          continue;
        for(size_t d=0; d<d_max; d++){
          double coef = (*blocks)[b].getRCoef(d,r);
          if(fabs(coef) > 1e-5)
            expr += coef*Yvar[b][t][d];
        }
      }
      if(cType == 'L'){
        model.addConstraint(expr <= probModel->getLimit(r, t)).setName("RL%d_%d",r,t);
      }else if(cType == 'G'){
        model.addConstraint(expr >= probModel->getLimit(r, t)).setName("RG%d_%d",r,t);
      }else{
        std::cerr << "ERROR: resource constraint type " << cType
          << " not implemented - IGNORED\n";
      }
    }
  }

  // objective function
  // TODO: check objective function, might have to go back and change all the
  //       variables to have their proper objective values


  if(lpFile != ""){
    model.writeLP(lpFile.c_str());
    std::cout << "Wrote " << lpFile << std::endl;
  }

  qol::Status status;
  bool timeout = true;

  while (timeout){
    try {
      timeout = false;
      status = solveRelaxed ? model.solveRelaxed() : model.solveExact();
    }
    catch (qol::Exception & ex) {
      std::cout << "Unable to find good solution in specified time, resolving!" << std::endl;
      timeout = true;
    }
  }

  std::cout << boost::format("Completed in %.2f sec CPU / %.2f sec wall. Objective = %f\n"
		       ) % timer.elapsedSeconds() % timer.elapsedWallTime() % (-model.getObjective());

  obj = -model.getObjective();

  std::string solnFile = probModel->getName()+"_window"+".sol";

  std::ofstream fp_out;
  fp_out.open(solnFile.c_str()); //, std::ios::app);

  size_t block_count = 0;
  std::vector<size_t> d_count(d_max,0);

  if(fp_out.is_open()){
    std::cout << "Writing solution to " << solnFile << std::endl;
    fp_out<< "# Status          "<<status<< std::endl;
    fp_out<< "# CPU time        "<<timer.elapsedSeconds()
      <<" sec.  Wall:"<< timer.elapsedWallTime() <<std::endl;
    fp_out<< "# Objective Value "<< obj << std::endl;
  }

  std::vector<int> update_count(d_max,0);

  //update Y_sol
  for (size_t b=0; b<nB; ++b){
    double b_sum = 0.0;
    for (size_t t =0; t<sh.WINDOW_SIZE; ++t){
      for (size_t d=0; d<d_max; ++d){
        if (model.getPrimal(Yvar[b][t][d]) != Y_sol[b][t+mod_tmin][d]){
          update_count[d]++;
        }

        if (model.getPrimal(Yvar[b][t][d]) > 1e-5){
          Y_sol[b][t+mod_tmin][d] = model.getPrimal(Yvar[b][t][d]);
        }
        else
          Y_sol[b][t+mod_tmin][d] = 0.0;

        b_sum += Y_sol[b][t+mod_tmin][d];
      }
    }
    if(b_sum - 1.0 > 0.00001){
      std::cout << std::setprecision(20) << "\n\t\tWarning 1:: Block: " << b
	         << " is used too often, sum: " << b_sum << ", violaton at: " << std::endl;
	    for(int t =0; t<t_max; t++){
	      for(int d =0; d<d_max; d++){
	        if(Y_sol[b][t][d] > 0.0){
            std::cout << "\t\t\tTime: " << t << ", dest: " << d << std::endl;
	        }
	      }
	    }
    }
  }

  std::cout << std::endl;
  for (size_t d=0;d<d_max;++d)
    std::cout << update_count[d] << " from dest " << d << " updated" << std::endl;

  std::cout << std::endl;
  delete mipPtrW;
  return obj;

}

#include "GRASPSolver.h"

void GRASPSolver::solveWindow(Sol_Int &sol, int t_0){
	
	std::vector<Block> * blocks=probModel->getBlock();

	std::vector<int> in_window(nB, -1);
	std::vector<int> window_map;
	std::vector<int> fixed;

	std::cout << "\nfinding blocks in window " << t_0 << " to " << t_0 + sh.WINDOW_SIZE-1 << ".. "<< std::flush;
	for (int b = 0; b < nB; ++b){
		if (sol.x[b] >= t_0 && sol.x[b] < t_0 + sh.WINDOW_SIZE){
			in_window[b] = window_map.size();
			window_map.push_back(b);

			// check if it is on the boundary with successors in the solution, if so, block must be fixed
			bool fix = false;
			std::vector<Node*> succ_nodes;
			probModel->graph.getNode(b)->getConnected(1,succ_nodes);
			// std::cout << "block " << b << " has " << succ_nodes.size() << " successors" <<std::endl;
			for (int succ_node_idx = 0; succ_nodes.size(); ++succ_node_idx){
				if (!include[succ_nodes[succ_node_idx]->getID()]){
					continue;
				}
				if (!(sol.x[succ_nodes[succ_node_idx]->getID()] >= t_0 && sol.x[succ_nodes[succ_node_idx]->getID()] < t_0 + sh.WINDOW_SIZE) && sol.x[succ_nodes[succ_node_idx]->getID()] < nT){
					fix = true;
					break;
				}
			}
			fixed.push_back(fix);
		}
	}
	std::cout << "done!\n" << std::endl;

	// set up vector for variables
	std::vector<std::vector<qol::Variable> > x = std::vector<std::vector< qol::Variable> >(window_map.size(), std::vector<qol::Variable> (sh.WINDOW_SIZE));

	qol::MIPSolver *mipPtr=0; // pointer to qol MIP solver


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
    // if (rel_gap)
      // param.setParamVal(qol::RELGAP,0.1);
    // else

    if (sh.MIP_GAP > 0)
      param.setParamVal(qol::RELGAP,sh.MIP_GAP);

    if (sh.WINDOW_SEARCH_TIME > 0)
      param.setParamVal(qol::TIMELIMIT,sh.WINDOW_SEARCH_TIME);



    // std::vector<bool> mined (nB,false);

    // determine solver type
    //if(solverType == GUROBI_T)
    //  mipPtr = new qol::GurobiFormulation();
    //else if(solverType == CPLEX_T)
    // if(solverType == CPLEX_T)
	  mipPtr = new qol::CplexFormulation();

    CPXENVptr env=(dynamic_cast<qol::CplexFormulation *>(mipPtr))->env;

    CPXsetintparam(env,CPX_PARAM_PREIND,sh.PRESOLVE);
    if (sh.HEURISTIC_SEARCH){
      CPXsetdblparam(env,CPX_PARAM_CUTSFACTOR,1.0);
      CPXsetintparam(env,CPX_PARAM_HEURFREQ,1);
      CPXsetintparam(env,CPX_PARAM_RINSHEUR,1);
      CPXsetintparam(env,CPX_PARAM_POLISHAFTERINTSOL,1);
    }    
    
    // create MIPSolver object
    qol::MIPSolver &mip=*mipPtr;
    mip.setParameters(param);

    //std::vector<std::vector<qol::Variable> > y(nB);

    std::cout << "Problem has " << nB << " blocks" << std::endl;

    initWindowModel(mip,x,in_window,window_map,sol,t_0, fixed);

    bool solveRelaxed = false;

    qol::CpuTimer timer;
    qol::Status status = solveRelaxed ? mip.solveRelaxed() : mip.solveExact();
    std::cout << boost::format("Completed in %.2f sec CPU / %.2f sec wall. Objective = %f\n"
			       ) % timer.elapsedSeconds() % timer.elapsedWallTime() % (-mip.getObjective());


    for (int b = 0;b<window_map.size();++b){
      for (int t = sh.WINDOW_SIZE-1; t >= 0; --t){
        if (mip.getPrimal(x[b][t]) > 1e-5){
          sol.x[window_map[b]] = t+t_0;
        }
      }
    }

    // sol.obj = 0.0;

    // for (int b = 0; b < nB; ++b){
    // 	if (sol.x[b] < nT){
    //     	sol.obj += (*blocks)[b].getProfit(0) / pow(1+rate,sol.x[b]);
    // 	}
    // }

    // std::cout << "\nTrue objective: " << sol.obj << std::endl;

    // std::cout << "\nChecking solution from solver...\n" << std::endl;
    //
    // bool test_error = verify((*probModel), sol);
    //
    // std::cout << "Solution found by solver was ";
    // if (!test_error)
    //   std::cout << "feasible!\n" << std::endl;
    // else{
    //   std::cout << "infeasible!\n" << std::endl;
    //   std::cin.get();
    //   throw qol::Exception("Infeasible solution!");
    // }

  //   std::ofstream fp_out;
  //   fp_out.open(solnFile.c_str()); //, std::ios::app);
  //   if(fp_out.is_open()){
  //     std::cout << "Writing solution to " << solnFile << std::endl;
  //     fp_out<< "# Status          "<<status<< std::endl;
  //     fp_out<< "# CPU time        "<<timer.elapsedSeconds()<<" sec.  Wall:"<< timer.elapsedWallTime() <<std::endl;
  //     fp_out<< "# Objective Value "<<-mip.getObjective() << std::endl;
  //     fp_out<< "# Block destination time yval"<<std::endl;
  //     /* solution file format blocks -> destinations -> time non-zero y value*/
  //     for(int b=0; b<nB; b++){
	// for(int t=0;t<t_max;++t)
	//   for(int d=0;d<d_max;++d){
	//     if(mip.getPrimal(y[b][d][t]) > 1e-5)
	//       fp_out << b << " " << d << " " << t << " " << mip.getPrimal(y[b][d][t]) << std::endl;
	//   }
  //     }
  //     fp_out.close();
    // }
   } // end try statement
  catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
  catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

  delete mipPtr;

}

void GRASPSolver::initWindowModel(qol::MIPSolver &mip, std::vector<std::vector<qol::Variable> > &x, const std::vector<int> &in_window, const std::vector<int> &window_map,Sol_Int &sol, int t_0,const std::vector<int> &fixed){

	std::vector<std::vector<int> > x_map (window_map.size(), std::vector<int> (sh.WINDOW_SIZE, -1));
	std::vector<Block> * blocks=probModel->getBlock();

	std::cout << "setting up decision variables.. " << std::flush;
	int var_count = 0;
    // decision variables
    for(int b_idx=0; b_idx<window_map.size(); b_idx++){

      double base_profit=(*blocks)[window_map[b_idx]].getProfit(0);
      for (int t=0; t<sh.WINDOW_SIZE; ++t){
        // compute profit for x[b][t] -> profit for block minus profit for block in next period
        double adjusted_profit = 0;
        if (t+t_0 < nT-1){
          adjusted_profit = base_profit/pow(1+rate,t+t_0) - base_profit/pow(1+rate,t+t_0+1);
        }
        else{
          adjusted_profit = base_profit/pow(1+rate,nT-1);
        }
        x[b_idx][t]=mip.addVar(0,1,-adjusted_profit,qol::Variable::BINARY, // LB,UB,cost
            boost::str(boost::format("x%03d_%03d")%b_idx%t));
        x_map[b_idx][t] = var_count++;
      }
    }

    std::vector<double> start_vec(window_map.size()*sh.WINDOW_SIZE,0.0);

    for (int b = 0; b<window_map.size(); ++b){
      for (int t = 0; t<sh.WINDOW_SIZE; ++t){
        if (sol.x[window_map[b]] <= t+t_0){
		    // mip.setPrimalStart(x[b][t],1);
	        start_vec[x_map[b][t]] = 1;
        }
      }
    }


    std::cout << "done!" << std::endl << "setting precedent constraints.. " <<std::flush;

    // type of constraints
    /*Precedence (7) (but using different definition of x */
    for(int a=0; a<window_map.size(); a++){
      std::vector<int> * pred = (*blocks)[window_map[a]].getPreds();
      int n = (*blocks)[window_map[a]].getNumPred();
      for(int p=0; p<n; p++){
        int b = (*pred)[p];
        if (in_window[b] > -1){
	        for (int t=0; t<sh.WINDOW_SIZE; ++t){
	          mip.addConstraint(x[in_window[b]][t] >= x[a][t]).setName("pre%d_%d_%d",b,a,t);
	        }
    	}
      }
    }

    std::cout << "done!" << std::endl << "setting remains done constraints.. " <<std::flush;

    /* Block remains done (9) */
    for(int b=0; b<window_map.size(); b++){
      for(int t=1; t<sh.WINDOW_SIZE; t++){
        mip.addConstraint(x[b][t-1] <= x[b][t]).setName("xt%d_%d",b,t);
      }
    }
    std::cout << "done!" << std::endl << "fixing boundary blocks.. "<<std::flush;

    for(size_t b=0; b<window_map.size(); b++){
		if (fixed[b]){
			mip.setVarLB(x[b][sh.WINDOW_SIZE-1],1);
			start_vec[x_map[b][sh.WINDOW_SIZE-1]] = 1;
		    // qol::Expression xbt_minus1 = x[b][sh.WINDOW_SIZE-1];
		    // mip.addConstraint(xbt_minus1 == 1).setName("lxw_%d",b);
	  	}
	}
	std::cout << "done!" << std::endl;
    

    // last window constraint to ensure every block completes
    if (t_0 == (nT-1)-sh.WINDOW_SIZE){
    	std::cout << "setting last window bounds.. " << std::flush;
		for(size_t b=0; b<window_map.size(); b++){
		//     qol::Expression xbt_minus1 = x[b][sh.WINDOW_SIZE-1];
		//     mip.addConstraint( xbt_minus1 == 1).setName("lxw_%d",b);
			mip.setVarLB(x[b][sh.WINDOW_SIZE-1],1);
			start_vec[x_map[b][sh.WINDOW_SIZE-1]] = 1;
		}
		std::cout << "done!" << std::endl;
    }

    mip.setPrimalStart(start_vec);
    
    std::cout << "setting resource constraints.. " << std::flush;

    /* Resource constraints (10) */
    for(int r=0; r<nR; r++){
      for (int t=0; t<sh.WINDOW_SIZE; ++t){
        char cType = probModel->getResConstrType(r, t+t_0);
        qol::Expression expr;
        for(int b=0; b<window_map.size(); b++){
          double coef = (*blocks)[window_map[b]].getRCoef(0,r);
          if(fabs(coef) > 1e-5){
            if (t > 0){
              expr += coef*(x[b][t]-x[b][t-1]);
            }
            else{
              expr += coef*x[b][t];
            }
          }
        }
        if(cType == 'L'){
          mip.addConstraint(expr <= probModel->getLimit(r,t)).setName("R%d_%d",r,t);
        }else if(cType == 'G'){
          mip.addConstraint(expr >= probModel->getLimit(r,t)).setName("R%d_%d",r,t);
        }else{
          std::cerr << "ERROR: resource constraint type " << cType
            << " not implemented - IGNORED\n";
        }
      }
    }
    std::cout << "done!\n" << std::endl;
}
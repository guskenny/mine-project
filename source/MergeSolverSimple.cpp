#include "MergeSolverSimple.h"

MergeSolverSimple::MergeSolverSimple(const SettingsHandler sh, SinglePModel *probModel, Sol_Int &best_sol, std::vector<std::vector<int> > &groups, std::vector<int> &group_map, std::vector<int> &fixed) : probModel(probModel), best_sol(best_sol), sh(sh),  groups(groups), group_map(group_map), fixed(fixed) {

  nB = probModel->getNBlock();
  nG = groups.size();
  nR = probModel->getnResources();
  nT = probModel->getNPeriod();
  rate = probModel->getDiscountRate();
  rel_gap = true;

  // set up x[b][t] vector for variables
  x = std::vector<std::vector< qol::Variable> >(nB, std::vector<qol::Variable> (nT));
  y = std::vector<std::vector< qol::Variable> >(nB, std::vector<qol::Variable> (nT));

  // initMergeModel();

  std::cout << "\nMerge Solver initialised" << std::endl;
}



void MergeSolverSimple::solve(Sol_Int &sol){
  qol::CpuTimer timer;

  qol::MIPSolver *mipPtr=0; // pointer to qol MIP solver

  std::vector<Block> * blocks=probModel->getBlock();

  std::string lpFile="";// ,solnFile = probModel->getName()+".sol";

  double bestObj = -1e9;

  //  for (int iter = 0; iter < sh.NUM_ITER; ++iter){
  //    std::cout << "\n******* ITERATION " << iter << " of "
  //              << sh.NUM_ITER << " *******\n" << std::endl;
  double trueObj = 0;

  rel_gap = false;

   try{
    // set verbosity
    qol::Parameters param;
    param.setParamVal(qol::VERBOSITY,1);
    if (rel_gap)
      param.setParamVal(qol::RELGAP,0.1);
    // else
    //   param.setParamVal(qol::RELGAP,0.001);

    if (sh.WINDOW_SEARCH_TIME > 0)
      param.setParamVal(qol::TIMELIMIT,sh.WINDOW_SEARCH_TIME);



    // std::vector<bool> mined (nB,false);

    // determine solver type
    //if(solverType == GUROBI_T)
    //  mipPtr = new qol::GurobiFormulation();
    //else if(solverType == CPLEX_T)
    // if(solverType == CPLEX_T)
	  mipPtr = new qol::CplexFormulation();

    // create MIPSolver object
    qol::MIPSolver &mip=*mipPtr;
    mip.setParameters(param);

    //std::vector<std::vector<qol::Variable> > y(nB);

    std::cout << "Problem has " << nB << " blocks" << std::endl;


    // initMergeModel(mip);

    initMergeSimpleModel(mip);

    bool solveRelaxed = false;

    qol::CpuTimer timer;
    qol::Status status = solveRelaxed ? mip.solveRelaxed() : mip.solveExact();
    std::cout << boost::format("Completed in %.2f sec CPU / %.2f sec wall. Objective = %f\n"
			       ) % timer.elapsedSeconds() % timer.elapsedWallTime() % (-mip.getObjective());

    // test solution for feasibility
    sol.init(nB,1);
    sol.obj = -mip.getObjective();
    sol.nT = nT;

    // for (int b = 0;b<nB;++b){
    //   for (int t = 0; t<nT;++t){
    //     sol.x[b] = -1;
    //   }
    // }

    for (int b = 0;b<nB;++b){
      int min_period = nT;
      for (int t = 0; t<nT;++t){
        if (mip.getPrimal(x[b][t]) > 1e-5){
          min_period = std::min(min_period,t);
        }
      }
      sol.x[b] = min_period;
    }

    // for (int b = 0;b<nB;++b){
    //   for (int t = 0; t<nT;++t){
    //     if (mip.getPrimal(y[b][t]) > 1e-5){
    //       sol.x[b] = t;
    //     }
    //   }
    // }

    std::cout << "\nTrue objective: " << sol.obj << std::endl;

  //   std::cout << "\nChecking solution from solver...\n" << std::endl;
  //
  //   bool test_error = verify((*probModel), sol);
  //
  //   std::cout << "Solution found by solver was ";
  //   if (!test_error)
  //     std::cout << "feasible!\n" << std::endl;
  //   else{
  //     std::cout << "infeasible!\n" << std::endl;
  //     throw qol::Exception("Infeasible solution!");
  //   }
   } // end try statement
  catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
  catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

  delete mipPtr;
}

//TODO: revert back to only single variables
void MergeSolverSimple::initMergeSimpleModel(qol::MIPSolver &mip){

    std::vector<Block> * blocks=probModel->getBlock();

    std::cout << "setting up decision variables.. " << std::flush;
    // decision variables
    for(int b=0; b<nB; b++){
      double base_profit=(*blocks)[b].getProfit(0);
      for (int t=0; t<nT; ++t){
        // compute profit for x[b][t] -> profit for block minus profit for block in next period
        double adjusted_profit = 0;
        if (t < nT-1){
          adjusted_profit = base_profit/pow(1+rate,t) - base_profit/pow(1+rate,t+1);
        }
        else{
          adjusted_profit = base_profit/pow(1+rate,nT-1);
        }
        x[b][t]=mip.addVar(0,0,-adjusted_profit,qol::Variable::BINARY, // LB,UB,cost
            boost::str(boost::format("x%03d_%03d")%b%t));
      }
    }

    std::cout << "done!" << std::endl << "setting bounds.. " <<std::flush;
    // unfix all variables in current window that are to be included

    int unfixed = 0;
    int one_fix = 0;

    for (int b = 0; b<nB; ++b){
      for (int t = 0; t<nT; ++t){
        if (fixed[t*nB + b] != 0){
          mip.setVarUB(x[b][t],1);
          if (fixed[t*nB + b] > 0){
            mip.setVarLB(x[b][t],1);
            one_fix++;
          }
          else{
            unfixed++;
          }
        }
        if (best_sol.x[b] <= t){
          mip.setPrimalStart(x[b][t],1);
        }
      }
    }

    // for (int b = 0; b<nB; ++b){
    //   if (best_sol.x[b] < 0){
    //     continue;
    //   }
    //   for (int t = 0; t < nT; ++t){
    //     if (t > merged.fixed[b][0]){
    //       mip.setVarUB(x[b][t],1);
    //       unfixed++;
    //     }
    //     if (t > merged.fixed[b][1]){
    //       mip.setVarLB(x[b][t],1);
    //       one_fix++;
    //     }
    //     if (best_sol.x[b] <= t){
    //       mip.setPrimalStart(x[b][t],1);
    //     }
    //   }
    // }

    std::cout << "done!" << std::endl << unfixed << " variables un-fixed from zero"
              << std::endl << one_fix << " variables fixed to one" <<std::endl;
    std::cout << "percent reduction in variables: " << double(unfixed/(nB*nT)) << std::endl << std::endl;
    std::cout << "setting up constraints:" << std::endl << "precedence.. " << std::flush;

    int skipped = 0;

    // get successor list
    std::vector<std::vector<int> > succ(nB, std::vector<int>());
    for (int b=0; b<nB; ++b){
      std::vector<int> * pred = (*blocks)[b].getPreds();
      int n = (*blocks)[b].getNumPred();
      for(int p=0; p<n; p++){
        int b_prev = (*pred)[p];
        succ[b_prev].push_back(b);
      }
    }

    // put in precedence constraints using successors
    for (int b=0; b<nB; ++b){
      for (std::vector<int>::iterator it = succ[b].begin(); it != succ[b].end(); ++it){
        for (int t=0; t < nT; ++t){
          if (fixed[t*nB+b] != -1 || fixed[t*nB + (*it)] == 0){
            skipped++;
            continue;
          }
          mip.addConstraint(x[b][t] >= x[*it][t]).setName("pre%d_%d_%d",b,*it,t);
        }
      }
    }

    //
    // // type of constraints
    // /*Precedence (7) (but using different definition of x */
    // for(int a=0; a<nB; a++){
    //   std::vector<int> * pred = (*blocks)[a].getPreds();
    //   int n = (*blocks)[a].getNumPred();
    //   for(int p=0; p<n; p++){
    //     int b = (*pred)[p];
    //     for (int t=0; t<nT; ++t){
    //       mip.addConstraint(x[b][t] >= x[a][t]).setName("pre%d_%d_%d",b,a,t);
    //     }
    //   }
    // }
    std::cout << "done! (" << skipped << " precedence vars skipped)\n" << std::endl << "block remains done.. " <<std::flush;

    int sum_skipped = 0;
    skipped = 0;
    /* Block remains done (9) */
    for(int b=0; b<nB; b++){
      for(int t=1; t<nT; t++){
        if (fixed[(t-1)*nB+b] == 0 || fixed[t*nB+b] == 1){
          skipped++;
          continue;
        }
        mip.addConstraint( x[b][t-1] <= x[b][t]).setName("xt%d_%d",b,t);
      }
    }
    std::cout << "done! (" << skipped << " remain done vars skipped)\n" << std::endl << "resource.. " <<std::flush;
    sum_skipped += skipped;
    skipped = 0;
    /* Resource constraints (10) */
    for(int r=0; r<nR; r++){
      for (int t=0; t<nT; ++t){
        char cType = probModel->getResConstrType(r, t);
        qol::Expression expr;
        for(int b=0; b<nB; b++){
          if (fixed[t*nB + b] == 0){ // if fixed to 0, skip adding to resources
            skipped++;
            continue;
          }
          double coef = (*blocks)[b].getRCoef(0,r);
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

    std::cout << "done! (" << skipped << " resource vars skipped)\n" << std::endl;

    std::cout << "\ntotal constraints skipped: " << sum_skipped << std::endl << std::endl;

    // std::cout << "done!" << std::endl << "aggregated groups.. " <<std::flush;
    //
    // for (int g = 0; g < groups.size(); ++g){
    //   if (groups[g].size()==1){
    //     continue;
    //   }
    //   // make pairs of all variables in groups
    //   for(int i = 0; i < groups[g].size(); ++i){
    //     qol::Expression lhs = x[groups[g][i]%nB][groups[g][i]/nB];
    //     for (int j = i+1; j < groups[g].size()-1; ++j){
    //       qol::Expression rhs = x[groups[g][j]%nB][groups[g][j]/nB];
    //       mip.addConstraint(lhs == rhs).setName("G%d",g);
    //     }
    //   }
    // }
    //
    // // for (int g = 0; g < groups.size(); ++g){
    // //   // if one selected, then all must be selected
    // //   qol::Expression if_selected = x[groups[g][0]%nB][groups[g][0]/nB]*int(groups[g].size());
    // //   qol::Expression sumX;
    // //   for(int v = 0; v < groups[g].size(); ++v){
    // //     sumX += x[groups[g][v]%nB][groups[g][v]/nB];
    // //   }
    // //   mip.addConstraint(sumX == if_selected).setName("G%d",g);
    // // }
    // std::cout << "done!" << std::endl;
}

//TODO: revert back to only single variables
void MergeSolverSimple::initMergeModel(qol::MIPSolver &mip){

    std::vector<Block> * blocks=probModel->getBlock();

    std::cout << "setting up decision variables.. " << std::flush;
    // decision variables
    for(int b=0; b<nB; b++){
      double base_profit=(*blocks)[b].getProfit(0);
      for (int t=0; t<nT; ++t){
        // compute profit for x[b][t] -> profit for block minus profit for block in next period
        double adjusted_profit = 0;
        if (t < nT-1){
          adjusted_profit = base_profit/pow(1+rate,t) - base_profit/pow(1+rate,t+1);
        }
        else{
          adjusted_profit = base_profit/pow(1+rate,nT-1);
        }
        x[b][t]=mip.addVar(0,0,-adjusted_profit,qol::Variable::BINARY, // LB,UB,cost
            boost::str(boost::format("x%03d_%03d")%b%t));
      }
    }

    std::cout << "done!" << std::endl << "setting bounds.. " <<std::flush;
    // unfix all variables in current window that are to be included

    int unfixed = 0;
    int one_fix = 0;
    int zero_fix = nB*nT;

    std::vector<int> matches(3,0);

    // for (int b = 0; b<nB; ++b){
    //   for (int t = 0; t < nT; ++t){
    //     if (t >= merged.fixed[b][0]){
    //       mip.setVarUB(x[b][t],1);
    //       if (t > merged.fixed[b][1]){
    //         if(fixed[t*nB + b] == 1){
    //           matches[1]++;
    //         }
    //         else{
    //           std::cout << "x[" << b << "][" << t << "] not matching, is "<< fixed[t*nB + b] << " should be 1 - [" << merged.fixed[b][0] << ", " << merged.fixed[b][1] << "]" << std::endl;
    //         }
    //         mip.setVarLB(x[b][t],1);
    //         one_fix++;
    //         zero_fix--;
    //
    //       }
    //       else{
    //         if(fixed[t*nB + b] == -1){
    //           matches[3]++;
    //         }
    //         else{
    //           std::cout << "x[" << b << "][" << t << "] not matching, is "<< fixed[t*nB + b] << " should be uf - [" << merged.fixed[b][0] << ", " << merged.fixed[b][1] << "]" << std::endl;
    //         }
    //         unfixed++;
    //         zero_fix--;
    //       }
    //     }
    //     else{
    //       if(fixed[t*nB + b] == 0){
    //         matches[0]++;
    //       }
    //       else{
    //         std::cout << "x[" << b << "][" << t << "] not matching, is "<< fixed[t*nB + b] << " should be 0 - [" << merged.fixed[b][0] << ", " << merged.fixed[b][1] << "]" << std::endl;
    //       }
    //     }
    //     if (best_sol.x[b] <= t){
    //       mip.setPrimalStart(x[b][t],1);
    //     }
    //   }
    // }

    std::cout << "done!" << std::endl << zero_fix << " variables fixed to 0" << std::endl << one_fix << " variables fixed to one"
              << std::endl << unfixed << " variables un-fixed from zero"
              << std::endl;
    std::cout << "setting up constraints:" << std::endl << "precedence.. " << std::flush;

    // type of constraints
    /*Precedence (7) (but using different definition of x */
    for(int a=0; a<nB; a++){
      std::vector<int> * pred = (*blocks)[a].getPreds();
      int n = (*blocks)[a].getNumPred();
      for(int p=0; p<n; p++){
        int b = (*pred)[p];
        for (int t=0; t<nT; ++t){
          mip.addConstraint(x[b][t] >= x[a][t]).setName("pre%d_%d_%d",b,a,t);
        }
      }
    }
    std::cout << "done!" << std::endl << "block remains done.. " <<std::flush;

    /* Block remains done (9) */
    for(int b=0; b<nB; b++){
      for(int t=1; t<nT; t++){
        mip.addConstraint( x[b][t-1] <= x[b][t]).setName("xt%d_%d",b,t);
      }
    }
    std::cout << "done!" << std::endl << "resource.. " <<std::flush;

    /* Resource constraints (10) */
    for(int r=0; r<nR; r++){
      for (int t=0; t<nT; ++t){
        char cType = probModel->getResConstrType(r, t);
        qol::Expression expr;
        for(int b=0; b<nB; b++){
          double coef = (*blocks)[b].getRCoef(0,r);
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

    // std::cout << "done!" << std::endl << "aggregated groups.. " <<std::flush;
    //
    // for (int g = 0; g < groups.size(); ++g){
    //   if (groups[g].size()==1){
    //     continue;
    //   }
    //   // make pairs of all variables in groups
    //   for(int i = 0; i < groups[g].size(); ++i){
    //     qol::Expression lhs = x[groups[g][i]%nB][groups[g][i]/nB];
    //     for (int j = i+1; j < groups[g].size()-1; ++j){
    //       qol::Expression rhs = x[groups[g][j]%nB][groups[g][j]/nB];
    //       mip.addConstraint(lhs == rhs).setName("G%d",g);
    //     }
    //   }
    // }
    //
    // // for (int g = 0; g < groups.size(); ++g){
    // //   // if one selected, then all must be selected
    // //   qol::Expression if_selected = x[groups[g][0]%nB][groups[g][0]/nB]*int(groups[g].size());
    // //   qol::Expression sumX;
    // //   for(int v = 0; v < groups[g].size(); ++v){
    // //     sumX += x[groups[g][v]%nB][groups[g][v]/nB];
    // //   }
    // //   mip.addConstraint(sumX == if_selected).setName("G%d",g);
    // // }
    // std::cout << "done!" << std::endl;
}


//
// void MergeSolverSimple::initMergeModel(qol::MIPSolver &mip){
//
//     std::vector<Block> * blocks=probModel->getBlock();
//
//     std::cout << "setting up decision variables.. " << std::flush;
//     // decision variables
//     for(int b=0; b<nB; b++){
//       for (int t=0; t<nT; ++t){
//         double profit=(*blocks)[b].getProfit(0)/pow(1+rate,t);
//         // compute profit for x[b][t] -> profit for block minus profit for block in next period
//         // double adjusted_profit = 0;
//         // if (t < nT-1){
//         //   adjusted_profit = base_profit/pow(1+rate,t) - base_profit/pow(1+rate,t+1);
//         // }
//         // else{
//         //   adjusted_profit = base_profit/pow(1+rate,nT-1);
//         // }
//         x[b][t]=mip.addVar(0,1,0.0,qol::Variable::BINARY, // LB,UB,cost
//             boost::str(boost::format("x%03d_%03d")%b%t));
//         y[b][t]=mip.addVar(0,1,-profit,qol::Variable::BINARY, // LB,UB,cost
//             boost::str(boost::format("y%03d_%03d")%b%t));
//       }
//     }
//
//     std::cout << "done!" << std::endl << "setting bounds.. " <<std::flush;
//     // unfix all variables in current window that are to be included
//
//     for (int b = 0; b<nB; ++b){
//       for (int t = 0; t < nT; ++t){
//         // if (t > merged.fixed[b][0]){
//         //   mip.setVarUB(x[b][t],1);
//         //   mip.setVarUB(y[b][t],1);
//         // }
//         // if (t > merged.fixed[b][1]){
//         //   mip.setVarLB(x[b][t],1);
//         //   mip.setVarLB(y[b][t],1);
//         // }
//         if (best_sol.x[b] <= t){
//           mip.setPrimalStart(x[b][t],1);
//         }
//       }
//     }
//
//     std::cout << "done!" << std::endl;
//     std::cout << "setting up constraints:" << std::endl << "precedence.. " << std::flush;
//
//     // type of constraints
//     /*Precedence (7) (but using different definition of x */
//     for(int a=0; a<nB; a++){
//       std::vector<int> * pred = (*blocks)[a].getPreds();
//       int n = (*blocks)[a].getNumPred();
//       for(int p=0; p<n; p++){
//         int b = (*pred)[p];
//         for (int t=0; t<nT; ++t){
//           mip.addConstraint(x[b][t] >= x[a][t]).setName("pre%d_%d_%d",b,a,t);
//         }
//       }
//     }
//     std::cout << "done!" << std::endl << "block remains done.. " <<std::flush;
//
//     /* SumDest (8) */
//     for(int b=0; b<nB; b++){
//       for(int t=0; t<nT; t++){
//         qol::Expression lhs = y[b][t];
//         qol::Expression rhs = x[b][t];
//         if(t > 0){
//           rhs -= x[b][t-1];
//         }
//         mip.addConstraint( lhs == rhs ).setName("dest%d_%d",b,t);
//       }
//     }
//
//     /* Block remains done (9) */
//     for(int b=0; b<nB; b++){
//       for(int t=1; t<nT; t++){
//         mip.addConstraint(x[b][t-1] <= x[b][t]).setName("xt%d_%d",b,t);
//       }
//     }
//     std::cout << "done!" << std::endl << "resource.. " <<std::flush;
//
//     /* Resource constraints (10) */
//     for(int r=0; r<nR; r++){
//       for (int t=0; t<nT; ++t){
//         char cType = probModel->getResConstrType(r, t);
//         qol::Expression expr;
//         for(int b=0; b<nB; b++){
//           double coef = (*blocks)[b].getRCoef(0,r);
//           if(fabs(coef) > 1e-5){
//             expr += coef*y[b][t];
//           }
//         }
//         if(cType == 'L'){
//           mip.addConstraint(expr <= probModel->getLimit(r,t)).setName("R%d_%d",r,t);
//         }else if(cType == 'G'){
//           mip.addConstraint(expr >= probModel->getLimit(r,t)).setName("R%d_%d",r,t);
//         }else{
//           std::cerr << "ERROR: resource constraint type " << cType
//             << " not implemented - IGNORED\n";
//         }
//       }
//     }
//     std::cout << "done!" << std::endl << "aggregated groups.. " <<std::flush;
//
//     // for (int g = 0; g < groups.size(); ++g){
//     //   // if one selected, then all must be selected
//     //   qol::Expression if_selected = x[groups[g][0]%nB][groups[g][0]/nB]*int(groups[g].size());
//     //   qol::Expression sumX;
//     //   for(int v = 0; v < groups[g].size(); ++v){
//     //     sumX += x[groups[g][v]%nB][groups[g][v]/nB];
//     //   }
//     //   mip.addConstraint(sumX == if_selected).setName("G%d",g);
//     // }
//     std::cout << "done!" << std::endl;
// }

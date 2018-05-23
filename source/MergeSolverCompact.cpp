#include "MergeSolverCompact.h"

MergeSolverCompact::MergeSolverCompact(const SettingsHandler sh, SinglePModel *probModel, Sol_Int &best_sol, std::vector<std::vector<int> > &groups, std::vector<int> &group_map, std::vector<int> &fixed, const std::vector<int> &include, std::ofstream *red_data) : probModel(probModel), best_sol(best_sol), sh(sh),  groups(groups), group_map(group_map), fixed(fixed), include(include), red_data(red_data) {

  nB = probModel->getNBlock();
  nG = groups.size();
  nR = probModel->getnResources();
  nT = probModel->getNPeriod();
  rate = probModel->getDiscountRate();
  rel_gap = true;

  x_map = std::vector<std::vector<int> > (nB, std::vector<int> (nT, -1));

  // int count = 0;
  // for (int b = 0; b < nB; ++b){
  //   for (int t = 0; t < nT; ++t){
  //     if (fixed[t*nB + b] != 0){
  //       x_map[b][t] = count++;
  //     }
  //   }
  // }
  //
  // std::cout << std::endl << count << " non-zero variables found" << std::endl;
  //
  // // set up x vector for variables
  // x = std::vector<qol::Variable> (count);

  // initMergeModel();

  std::cout << "\nMerge Solver initialised" << std::endl;
}


void MergeSolverCompact::solve(Sol_Int &sol){
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
    // if (rel_gap)
      // param.setParamVal(qol::RELGAP,0.1);
    // else
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

    // create MIPSolver object
    qol::MIPSolver &mip=*mipPtr;
    mip.setParameters(param);

    //std::vector<std::vector<qol::Variable> > y(nB);

    std::cout << "Problem has " << nB << " blocks" << std::endl;


    // initMergeModel(mip);

    if (sh.GROUP_MERGE){
      initMergeGroupModel(mip);
    }
    else{
      initMergeSimpleModel(mip);
    }

    bool solveRelaxed = false;

    qol::CpuTimer timer;
    qol::Status status = solveRelaxed ? mip.solveRelaxed() : mip.solveExact();
    std::cout << boost::format("Completed in %.2f sec CPU / %.2f sec wall. Objective = %f\n"
			       ) % timer.elapsedSeconds() % timer.elapsedWallTime() % (-mip.getObjective());


   if (sh.GROUP_MERGE){
     extractGroupSol(mip, sol);
   }
   else{
     extractSimpleSol(mip, sol);
   }

    // for (int b = 0;b<nB;++b){
    //   for (int t = 0; t<nT;++t){
    //     if (mip.getPrimal(y[b][t]) > 1e-5){
    //       sol.x[b] = t;
    //     }
    //   }
    // }

    std::cout << "\nTrue objective: " << sol.obj << std::endl;

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

void MergeSolverCompact::extractSimpleSol(qol::MIPSolver &mip, Sol_Int &sol){
  // test solution for feasibility
  sol.init(nB,nT);
  sol.obj = -mip.getObjective();
  sol.nT = nT;

  // for (int b = 0;b<nB;++b){
  //   for (int t = 0; t<nT;++t){
  //     sol.x[b] = -1;
  //   }
  // }

  for (int b = 0;b<nB;++b){
    if (!include[b]){
      continue;
    }
    int min_period = nT;
    for (int t = 0; t<nT;++t){
      if (x_map[b][t] < 0){
        continue;
      }
      if (mip.getPrimal(x[x_map[b][t]]) > 1e-5){
        min_period = std::min(min_period,t);
      }
    }
    sol.x[b] = min_period;
  }

  // // test solution for feasibility
  // sol.init(nB,nT);
  // sol.obj = -mip.getObjective();
  //
  // // for (int b = 0;b<nB;++b){
  // //   for (int t = 0; t<nT;++t){
  // //     sol.x[b] = -1;
  // //   }
  // // }
  //
  // for (int b = 0;b<nB;++b){
  //   // int min_period = nT;
  //   for (int t = 0; t<nT;++t){
  //     if (x_map[b][t] < 0){
  //       continue;
  //     }
  //     if (mip.getPrimal(x[x_map[b][t]]) > 1e-5){
  //       sol.x[b] = std::min(sol.x[b],t);
  //     }
  //   }
  //   // sol.x[b] = min_period;
  // }
}

void MergeSolverCompact::extractGroupSol(qol::MIPSolver &mip, Sol_Int &sol){
  // test solution for feasibility
  sol.init(nB,nT);
  sol.obj = -mip.getObjective();

  int groups_selected = 0;

  for (int b = 0; b < nB; ++b){
    if (!include[b]){
      continue;
    }
    for (int t = 0; t < nT; ++t){
      if (mip.getPrimal(x[group_map[nB*t + b]])){
        sol.x[b] = std::min(sol.x[b],t);
      }
    }
  }



  //
  // for (int group = 0; group < nG; ++group){
  //   if (mip.getPrimal(x[group])){
  //     groups_selected++;
  //     for (int member_idx = 0; member_idx < groups[group].size(); ++member_idx){
  //       int member = groups[group][member_idx];
  //       int block_idx = member % nB;
  //       int period = member / nB;
  //       sol.x[block_idx] = std::min(sol.x[block_idx],period);
  //     }
  //   }
  // }
  std::cout << groups_selected << " groups selected of " << nG << " total groups" << std::endl;
}

void MergeSolverCompact::initMergeGroupModel(qol::MIPSolver &mip){
  std::cout << "*** creating aggregate merge MIP model ***" << std::endl;

  std::vector<Block> * blocks=probModel->getBlock();

  // for debugging
  int num_prec_constr = 0;
  int num_res_constr = 0;
  int num_block_constr = 0;

  // create adjacency matrix for precedence constraints so only one constraint is made per group pair
  std::vector<std::vector<int> > group_matrix(nG, std::vector<int> (nG,0));

  // create vector for cumulative group resource usages:
  // vector size is |nG| x |nT| x |nR|
  std::vector<std::vector<std::vector<long double> > > group_rCoef(nG, std::vector<std::vector<long double> >(nT, std::vector<long double>(nR, 0.0)));

  // create vector for group profits
  std::vector<long double> group_profits(groups.size());

  // create variables for each group, need to do this before so precedence
  // constraints can be created on the fly
  for (int group = 0; group < nG; ++group){
    x.push_back(mip.addVar(0,1,0,qol::Variable::BINARY,(boost::format("group_%s") %group).str()));
  }

  // iterate over all groups
  for (int group = 0; group < nG; ++group){
    // value to hold cumulative profit for group
    long double group_profit = 0;

    // iterate over each member (block/time pair) in the group
    for (int member = 0; member < groups[group].size(); ++member){
      // get current block id and block object
      int curr_member = groups[group][member];
      int block_idx = curr_member % nB; // gets block id
      int curr_period = curr_member/nB; // gets current period

      // add predecessor constraints
      std::vector<int> pred = (*blocks)[block_idx].getPred();
      for (int p = 0; p < pred.size(); ++p){
        // index of previous block
        int prev_block = pred[p];
        // index of block/time pair
        int prev_member = curr_period * nB + prev_block;
        int prev_group = group_map[prev_member];
        // add precedence constraint if one doesn't exist
        if (group != prev_group){ // no self arcs
          if (!group_matrix[prev_group][group]){
            mip.addConstraint(x[group] <= x[prev_group]).setName("pred_%d_%d",group,prev_group);
            num_prec_constr++;
            group_matrix[prev_group][group] = 1; // update matrix to avoid duplicates
          }
        }
      }

      // add arc to same block in previous period
      if (curr_period > 0){
        int prev_member = curr_member - nB;
        int prev_group = group_map[prev_member];
        // add constraint if one doesn't exist
        if (group != prev_group){ // no self arcs
          // reversed because arrows are pointing back in time
          if (!group_matrix[group][prev_group]){
            mip.addConstraint(x[prev_group] <= x[group]).setName("block_%d_%d",prev_group,group);
            num_block_constr++;
            group_matrix[group][prev_group] = 1; // update matrix
          }
        }
      }

      // calculate resource usage coefficient for member:
      // if a block occurs in the group then the set of time points for that block
      // within the group will be the interval [s, s+t, ..., t]
      // if group is selected, then block is mined at time s, so the resource coefficients
      // for that block will be added to the vector corresponding to period s
      // in order to cancel out subsequent groups, the resource coefficients for the block
      // must be subtracted from the vector corresponding to period t+1, this is because
      // in order for the solution to be feasible, the next group must start at t+1

      // calculate profit for each group:
      // similar to resource usage, add the profit block for the first time the block
      // appears in the group, and subtract the profit from the first period after last
      // period that block appears in the group so that the next group's profit will be
      // cancelled out

      // check if period member is the first time the block appears in the group
      if (curr_period == 0 || group != group_map[curr_member - nB]){
        // add current block profit for current period to current group
        long double block_profit=(*blocks)[block_idx].getProfit(0);
        group_profit += block_profit/pow(1+rate,curr_period);
        // add resources to resources vector for group in current period
        for (int r = 0; r < nR; ++r){
          group_rCoef[group][curr_period][r] += (*blocks)[block_idx].getRCoef(0, r);

        }
      }

      // check if current period is the last time the block appears in the group
      if (curr_period < nT-1){
        if (group != group_map[curr_member + nB]){
          // subtract current block profit for next period from current group
          long double block_profit=(*blocks)[block_idx].getProfit(0);
          group_profit -= block_profit/pow(1+rate,curr_period+1);
          for (int r = 0; r < nR; ++r){
            // subtract resources from resources vector in next period
            group_rCoef[group][curr_period+1][r] -= (*blocks)[block_idx].getRCoef(0, r);
          }
        }
      }
    }

    // set group cumulative profit
    mip.setObjCoeff(x[group], (double)-group_profit);
    group_profits[group] = group_profit;
  }

  // add resource constraints:
  // for each resource and period, the sum of all the coefficients of all
  // the groups must be less than the resource limit for that resource and period
  for (int r=0; r < nR; ++r){
    for (int t=0; t < nT; ++t){
      qol::Expression expr;
      for (int group=0; group < nG; ++group){
        double coef = (double)group_rCoef[group][t][r];
        if (fabs(coef) > 1e-5){
          expr += coef*x[group];
        }
      }
      mip.addConstraint(expr <= probModel->getLimit(r,t)).setName("R%d_%d",r,t);
      num_res_constr++;
    }
  }

  // set MIP start solution
  std::vector<double> start_vec(x.size(),0.0);

  if (sh.FIX_BEST_GROUP){
    int groups_fixed = 0;
    std::vector<int> stack;
    std::priority_queue<std::pair<double, int>> q;
    for (int g = 0; g < groups.size(); ++g) {
      q.push(std::pair<double, int>(group_profits[g]/groups[g].size(), g));
    }
    for (int i = 0; i < sh.FIX_BEST_GROUP; ++i) {
      int ki = q.top().second;
      double ks = q.top().first;
      std::cout << "fixing group " << ki << " group, " << ks << " = "
                << group_profits[ki] << " / " << groups[ki].size()
                << " variables" << std::endl;
      stack.push_back(ki);
      q.pop();
    }

    // set all groups before (in space) and after (in time) of best group to 1
    std::vector<int> visited(groups.size(), 0);
    while(!stack.empty()){
      int g = stack.back();
      groups_fixed++;
      stack.pop_back();
      if (visited[g]){
        continue;
      }
      visited[g] = 1;
      mip.setVarLB(x[g],1);
      start_vec[g] = 1;
      // iterate over the column associated with current group to get predecessors
      for (int i = 0; i < groups.size(); ++i){
        if (group_matrix[i][g] && !visited[i]){
          stack.push_back(i);
        }
      }
    }
    std::cout << groups_fixed << " groups fixed in total" << std::endl;
  }

  for (int b = 0; b<nB; ++b){
    for (int t = 0; t<nT; ++t){
      if (best_sol.x[b] <= t){
        // mip.setPrimalStart(x[x_map[b][t]],1);
        start_vec[group_map[nB * t + b]] = 1;
      }
    }
  }

  mip.setPrimalStart(start_vec);

  // fix all zero and all 1 group
  for (int i = 0; i < 2; ++i){
    for (int member = 0; member < fixed.size(); ++member){
      if (fixed[member] == i){
        mip.setVarUB(x[group_map[member]],i);
        mip.setVarLB(x[group_map[member]],i);
        break;
      }
    }
  }


  // for debugging
  std::cout << nG << " variables for groups added" << std::endl
            << num_prec_constr << " precedence constraints added" << std::endl
            << num_block_constr << " block persistence constraints added" << std::endl
            << num_res_constr << " resource constraints added" << std::endl;
}


//TODO: revert back to only single variables
void MergeSolverCompact::initMergeSimpleModel(qol::MIPSolver &mip){

    std::vector<Block> * blocks=probModel->getBlock();

    int one_fix = 0;
    int skipped = 0;
    int unfixed = 0;
    int var_count = 0;

    std::cout << "setting up decision variables.. " << std::flush;
    // decision variables
    for(int b=0; b<nB; b++){
      double base_profit=(*blocks)[b].getProfit(0);
      for (int t=0; t<nT; ++t){
        if (fixed[t*nB + b] == 0){
          skipped++;
          continue;
        }
        // compute profit for x[b][t] -> profit for block minus profit for block in next period
        double adjusted_profit = 0;
        if (t < nT-1){
          adjusted_profit = base_profit/pow(1+rate,t) - base_profit/pow(1+rate,t+1);
        }
        else{
          adjusted_profit = base_profit/pow(1+rate,nT-1);
        }
        if (fixed[t*nB + b] < 0){
          x.push_back(mip.addVar(0,1,-adjusted_profit,qol::Variable::BINARY, // LB,UB,cost
              boost::str(boost::format("x%03d_%03d")%b%t)));
          unfixed++;
        }
        else if (fixed[t*nB + b] > 0){
          x.push_back(mip.addVar(1,1,-adjusted_profit,qol::Variable::BINARY, // LB,UB,cost
              boost::str(boost::format("x%03d_%03d")%b%t)));
          one_fix++;
        }
        x_map[b][t] = var_count++;
      }
    }

    std::cout << "done!" << std::endl << "setting bounds.. " <<std::flush;
    // unfix all variables in current window that are to be included

    std::vector<double> start_vec(x.size(),0.0);

    for (int b = 0; b<nB; ++b){
      for (int t = 0; t<nT; ++t){
        if (best_sol.x[b] <= t){
          // mip.setPrimalStart(x[x_map[b][t]],1);
          start_vec[x_map[b][t]] = 1;
        }
      }
    }

    mip.setPrimalStart(start_vec);

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

    std::cout << "done!" << std::endl << std::endl << skipped << " variables skipped"
              << std::endl << one_fix << " variables fixed to one"
              << std::endl << unfixed << " variables unfixed"
              << std::endl << x.size() << " variables in MIP"
              << std::endl << nB*nT << " variables in problem" << std::endl;
    std::cout << "\nsetting up constraints:" << std::endl << "precedence.. " << std::flush;

    skipped = 0;
    int sum_skipped = 0;
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
          if (x_map[b][t] < 0 || x_map[*it][t] < 0 || fixed[t*nB + b] > 0){
            skipped++;
            continue;
          }
          mip.addConstraint(x[x_map[b][t]] >= x[x_map[*it][t]]).setName("pre%d_%d_%d",b,*it,t);
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
    std::cout << "done! (" << skipped << " precedence constraints skipped)\n" << std::endl << "block remains done.. " <<std::flush;

    sum_skipped += skipped;
    skipped = 0;
    /* Block remains done (9) */
    for(int b=0; b<nB; b++){
      for(int t=1; t<nT; t++){
        if (fixed[(t-1)*nB+b] == 0 || fixed[t*nB+b] == 1){
          skipped++;
          continue;
        }
        mip.addConstraint( x[x_map[b][t-1]] <= x[x_map[b][t]]).setName("xt%d_%d",b,t);
      }
    }
    std::cout << "done! (" << skipped << " remain done constraints skipped)\n" << std::endl << "resource.. " <<std::flush;

    sum_skipped += skipped;

    (*red_data) << unfixed << "," << sum_skipped << std::endl;

    skipped = 0;
    /* Resource constraints (10) */
    for(int r=0; r<nR; r++){
      for (int t=0; t<nT; ++t){
        qol::Expression expr;
        for(int b=0; b<nB; b++){
          if (fixed[t*nB + b] == 0){ // if fixed to 0, skip adding to resources
            skipped++;
            continue;
          }
          // if (t > 0){
          //   if (fixed[(t-1)*nB + b] == 0){ // if fixed to 0, skip adding to resources
          //     skipped++;
          //     continue;
          //   }
          // }
          double coef = (*blocks)[b].getRCoef(0,r);
          if(fabs(coef) > 1e-5){
            if (t > 0){
              if (fixed[(t-1)*nB + b] == 0){
                expr += coef*x[x_map[b][t]];
              }
              else{
                expr += coef*(x[x_map[b][t]]-x[x_map[b][t-1]]);
              }
            }
            else{
              expr += coef*x[x_map[b][t]];
            }
          }
        }
        // if(cType == 'L'){
          mip.addConstraint(expr <= probModel->getLimit(r,t)).setName("R%d_%d",r,t);
        // }else if(cType == 'G'){
        //   mip.addConstraint(expr >= probModel->getLimit(r,t)).setName("R%d_%d",r,t);
        // }else{
        //   std::cerr << "ERROR: resource constraint type " << cType
        //     << " not implemented - IGNORED\n";
        // }
      }
    }

    std::cout << "done! (" << skipped << " resource constraints skipped)\n" << std::endl;

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

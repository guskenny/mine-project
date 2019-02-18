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

void SinglePSolver::computeResUse(Sol_Int &sol){
  const int nB = probModel->getNBlock();
  const int t_max = probModel->getNPeriod();
  const int r_max = probModel->getnResources();
  const double rate = probModel->getDiscountRate();

  sol.res_use.clear();
  sol.res_use.resize(t_max, std::vector<long double>(r_max, 0.0));

  double mine_profit = 0;

  for (int b = 0; b < nB; ++b){
    if (sol.x[b] < t_max){
      const Block & block=probModel->getBlock(b);
      double block_profit =
      mine_profit += block.getProfit(0) / pow(1+rate,sol.x[b]);
      for (int r = 0; r < r_max; ++r){
        sol.res_use[sol.x[b]][r] += block.getRCoef(0,r);
      }
    }
  }
  sol.obj = mine_profit;
}

int SinglePSolver::solve(){
  switch(sh.MERGE_TYPE) {
    case 0 :{
       std::cout << "Running serialMergeSolve()\n" << std::endl;
       serialMergeSolve();
       break;
     }
    case 1 :{
       std::cout << "Running forkMergeSolve()\n" << std::endl;
       forkMergeSolve();
       break;
    }
    default :{
       std::cout << "Running randomMergeSolve()\n" << std::endl;
       // randomMergeSolve();
    }
  }
}

void SinglePSolver::computeUPIT(BranchNode_info &base_info, std::vector<int> &include){
  const int nB = probModel->getNBlock();
  const int t_max = probModel->getNPeriod();
  const int d_max = probModel->getnDestination();
  const int r_max = probModel->getnResources();
  std::vector<Block> * blocks=probModel->getBlock();

  //compute UPIT and remove blocks never mined
  UpitSolver upit((*probModel));
  upit.solve();
  upit.getClosure(include);
  int cnt=0;
  for(int b=0;b<nB;++b)
  if(!include[b]){
    base_info.time[b][0]=base_info.time[b][1]=t_max;
    if(++cnt < 10) std::cout << "\t\tBlock " << b << " is never mined\n";
    base_info.mine_block[b] = false;
  }
}

void SinglePSolver::getSeeds(BranchNode_info blank_info, std::vector<Sol_Int> &seeds){
  const int nB = probModel->getNBlock();
  const int t_max = probModel->getNPeriod();
  const int d_max = probModel->getnDestination();
  const int r_max = probModel->getnResources();
  std::vector<Block> * blocks=probModel->getBlock();

    Sol_Int best_sol(nB, d_max, r_max, t_max);
    double bestObj = -1e9;

    int n_periods = t_max;

    qol::CpuTimer run_timer;

    // establish vector for solutions
    std::vector<Sol_Int> init_seeds(sh.INIT_SEEDS);
    std::vector<double> it_times(sh.INIT_SEEDS);
    RandomSearch rs(sh);

    int best_seed = 0;

    double seed_start = run_timer.elapsedSeconds();

    #pragma omp parallel for
    for (int seed = 0; seed < sh.INIT_SEEDS; ++seed){
      std::cout << "\n******* SEED " << seed+1 << " of "
                << sh.INIT_SEEDS << " *******\n" << std::endl;

      Sol_Int sol_iter(nB, t_max);

      double before_timer = run_timer.elapsedSeconds();

      double iterObj = 0;

      if (sh.RANDOM_SEARCH){
        std::cout << "running randomSearch.. " << std::flush;
        n_periods = 0;
        iterObj = rs.randomSearch(blank_info, sol_iter, probModel); //TODO: CHANGE THIS PROPERLY!!
        std::cout << " done!" << std::endl;
      }
      else{
        ConeMiner cm = ConeMiner((*probModel),blank_info,sh);
        cm.solve(sol_iter);
        iterObj = sol_iter.obj;
      }

      // test final solution for feasibility
      try{
        // std::cout << "\nChecking solution from solver...\n" << std::endl;

        bool test_error = verify((*probModel), sol_iter);

        std::cout << "\nInitial solution is ";
        if (!test_error){
          std::cout << "feasible!\n" << std::endl;
        }
        else{
          std::cout << "\033[31;1minfeasible!\033[0m\n" << std::endl;
          std::cin.get();
          throw qol::Exception("Infeasible solution!");
          // sol = backup_sol;
        }
      } // end try statement
      catch (qol::Exception & ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
      }
      catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

      if (sh.TEST_SOL_MERGE){
        init_seeds[seed]=sol_iter;
      }

      double after_timer = run_timer.elapsedSeconds();

      it_times[seed] = after_timer-before_timer;



      if (iterObj > bestObj){
          std::cout << std::endl << "\033[32;1m***** NEW BEST OBJECTIVE FOUND! *****\033[0m\n" << std::endl;
          bestObj = iterObj;
          best_seed = seed;
          best_sol = sol_iter;
      }

    }//end of iteration loop

    std::cout << "\nTime to generate seeds: " << run_timer.elapsedSeconds() - seed_start << std::endl << std::endl;
    best_sol.obj = bestObj;

    // LocalSearch ls = LocalSearch(sh,probModel,include);

    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)

    // std::vector<Sol_Int> sols(merge_pop_size);
    // std::vector<Sol_Int> sols(sh.MERGE_POP_SIZE);

    computeResUse(best_sol);

    seeds.push_back(best_sol);

    std::uniform_int_distribution<int> uni(0,init_seeds.size()-1); // guaranteed unbiased

    // make sure all seeds are different
    std::vector<int> chosen_seeds;
    chosen_seeds.push_back(best_seed);

    std::cout << "adding seed " << best_seed << " to merge population" << std::endl;
    // for (int i = 0;i < sh.NUM_SEEDS-1; ++i){
    while (seeds.size() < sh.NUM_SEEDS){
      //add random seed in as well
      int rand_seed_idx = uni(rng); // index for open_idx
      if(std::find(chosen_seeds.begin(), chosen_seeds.end(), rand_seed_idx) != chosen_seeds.end()) {
        continue;
      }
      chosen_seeds.push_back(rand_seed_idx);
      Sol_Int temp_sol = init_seeds[rand_seed_idx];
      computeResUse(temp_sol);
      seeds.push_back(temp_sol);
      std::cout << "adding seed " << rand_seed_idx << " to merge population" << std::endl;
    }
}

void SinglePSolver::saveSols(const std::vector<Sol_Int> &sols, std::string path_name, int merge=0){
  const int nB = probModel->getNBlock();
  path_name = path_name + probModel->getName();
  for (int i = 0; i < sols.size(); ++i){
    std::ofstream solfile;
    std::string sol_idx = std::to_string(i);
    if (i < 10){
      sol_idx = "0" + std::to_string(i);
    }

    if (sols.size() == 1){
      sol_idx = "0" + std::to_string(merge);
    }

    std::string fname = path_name + "/" + probModel->getName() + "_" + sol_idx + ".sol";

    std::cout << "i: " << i <<  ", writing solution " << sol_idx << ", fname: " << fname << "\n";

    solfile.open(fname);
    for (int b = 0; b < nB; ++b){
      solfile << sols[i].x[b] << std::endl;
    }
    solfile.close();
  }
}

bool SinglePSolver::loadSols(std::vector<Sol_Int> &sols){
  const int nB = probModel->getNBlock();
  const int t_max = probModel->getNPeriod();
  const int r_max = probModel->getnResources();

  sols.clear();

  std::vector<double> obj_vals;

    for (int i = 0; i < sh.NUM_SEEDS; ++i){
    int cur_idx = sh.SOL_IDX + i;
    Sol_Int temp_sol(nB,0,r_max,t_max);

    std::string path_name = "./init_sols/" + probModel->getName();

    if (sh.RANDOM_SEARCH){
      path_name = "./rand_sols/" + probModel->getName();
    }
    
    // std::string sol_idx = std::to_string(i);

    // if (i < 10){
    //   sol_idx = "0" + std::to_string(i);
    // }

    std::string sol_idx = std::to_string(cur_idx);
    if (cur_idx < 10){
      sol_idx = "0" + std::to_string(cur_idx);
    }

    std::string fname = path_name + "/" + probModel->getName() + "_" + sol_idx + ".sol";
    std::ifstream solfile(fname);
    if (!solfile.is_open()){
      std::cout << "unable to find file: " << fname << "! Exiting...." << std::endl;
      exit(0);
    }
    std::vector<int> temp_x;
    std::string in_line;

    while(std::getline(solfile, in_line)){
      temp_x.push_back(stoi(in_line));
    }

    if (nB == temp_x.size()){
      temp_sol.x = temp_x;
    }
    else{
      std::cout << "error in solution file: " << fname << "! not the same number of entries!" << std::endl;
      return false;
    }

    std::cout << fname << " loaded!\n";

    computeResUse(temp_sol);
    sols.push_back(temp_sol);

    std::cout << "Adding solution " << cur_idx << " to seeds. Objective value: " << temp_sol.obj << std::endl;
    obj_vals.push_back(temp_sol.obj);

    solfile.close();
  }

  std::cout << sols.size() << " seeds loaded!" << std::endl;
  std::cout << "Best objective value: " << sols[0].obj << std::endl;
  std::cout << "Mean objective value of all seeds: " << my_math::calc_mean(obj_vals) << std::endl;
  std::cout << "Standard deviation: " << my_math::calc_std_dev(obj_vals) << std::endl;
  return true;
}

// set up MIP and solve
int SinglePSolver::forkMergeSolve(){

  qol::CpuTimer over_timer;

  const int nB = probModel->getNBlock();
  const int t_max = probModel->getNPeriod();
  const int d_max = probModel->getnDestination();
  const int r_max = probModel->getnResources();
  std::vector<Block> * blocks=probModel->getBlock();

  std::ofstream outfile;
  std::ofstream red_data;

  std::vector<std::ofstream> time_out(sh.NUM_SEEDS+1);

  std::vector<std::ofstream> seed_out((sh.NUM_SEEDS+1)*2);

  std::ofstream run_time;

  // std::ofstream sol_tracker;
  // sol_tracker.open("sol_tracker.csv");
  // sol_tracker << "SOLUTION TRACKER VALUES" << std::endl;
  // sol_tracker.close();
  //
  // std::ofstream sol_tracker_base;
  // sol_tracker_base.open("sol_tracker_base.csv");
  // sol_tracker_base << "SOLUTION TRACKER BASE VALUES" << std::endl;
  // sol_tracker_base.close();

  if (sh.RECORD_DATA){
  for (int i = 0; i < sh.NUM_SEEDS+1; ++i){
    time_out[i].open("./tracking_data/" + probModel->getName() + "/TIME_"+std::to_string(i)+".csv");
    time_out[i] << "TIME DATA - SEED " << i << std::endl << "CPU TIME,OBJ_VALUE" << std::endl;

    // output files for recording all solutions
    seed_out[i].open("./tracking_data/" + probModel->getName() + "/SEED_"+std::to_string(i)+".csv");
    seed_out[i] << "PARALLEL MERGE - SEED " << i << std::endl;
    seed_out[i] << "SA_ALPHA,SA_T_MIN,SA_ITER,FULL_RUNS,MERGE_POP_SIZE,NUM_ITER,ITER_INCR" << std::endl;
    seed_out[i] << sh.SA_ALPHA <<","<< sh.SA_T_MIN <<","<< sh.SA_ITER <<","<< sh.FULL_RUNS <<","<< sh.MERGE_POP_SIZE <<","<< sh.NUM_ITER << "," << sh.ITER_INCR << std::endl;

    // output files for recording best solutions
    seed_out[i+sh.NUM_SEEDS+1].open("./tracking_data/" + probModel->getName() + "/SEED_"+std::to_string(i)+"_BEST.csv");
    seed_out[i+sh.NUM_SEEDS+1] << "PARALLEL MERGE (BEST SOL)- SEED " << i << std::endl;
    seed_out[i+sh.NUM_SEEDS+1] << "SA_ALPHA,SA_T_MIN,SA_ITER,FULL_RUNS,MERGE_POP_SIZE,NUM_ITER,ITER_INCR" << std::endl;
    seed_out[i+sh.NUM_SEEDS+1] << sh.SA_ALPHA <<","<< sh.SA_T_MIN <<","<< sh.SA_ITER <<","<< sh.FULL_RUNS <<","<< sh.MERGE_POP_SIZE <<","<< sh.NUM_ITER << "," << sh.ITER_INCR << std::endl;
  }
  run_time.open ("./tracking_data/" + probModel->getName() + "/run_time.csv");
  run_time << "RUN_TIMES:" << std::endl;
  outfile.open ("./tracking_data/" + probModel->getName() + "/merge_obj.csv");
  red_data.open ("./tracking_data/" + probModel->getName() + "/reduce_data.csv");
  red_data << "VARIABLES,CONSTRAINTS" << std::endl;
  outfile << "PARALLEL MERGE" << std::endl;
  outfile << "SA_ALPHA,SA_T_MIN,SA_ITER,FULL_RUNS,MERGE_POP_SIZE,NUM_ITER,ITER_INCR" << std::endl;
  outfile << sh.SA_ALPHA <<","<< sh.SA_T_MIN <<","<< sh.SA_ITER <<","<< sh.FULL_RUNS <<","<< sh.MERGE_POP_SIZE <<","<< sh.NUM_ITER << "," << sh.ITER_INCR << std::endl;
  }

  BranchNode_info base_info(nB,t_max,d_max,r_max);
  std::vector<int> include(nB,0);

  computeUPIT(base_info, include);

  for (int full_run = 0; full_run < sh.NUM_TOTAL_RUNS; ++full_run){

    qol::CpuTimer run_timer;

    std::vector<Sol_Int> seeds;

    int merge_pop_size = sh.MERGE_POP_SIZE;

    int n_periods = t_max;

    std::string lpFile="";// ,solnFile = probModel->getName()+".sol";

    SolutionMerger sm(sh);
    LocalSearch ls = LocalSearch(sh,probModel,include);

    std::cout << probModel->getName() << std::endl;

    if (sh.LOAD_SEEDS){
      loadSols(seeds);
    }
    else{
      std::cout << "getting seeds.. "<<std::flush;
      getSeeds(base_info, seeds);
      std::cout << "done!" << std::endl;
    }

    Sol_Int best_sol = seeds[0];

    // Sol_Int merged_sol;

    // doMerge(sm, seeds, include, seeds[0] ,merged_sol,red_data);

    double bestObj = best_sol.obj;

    std::vector<Sol_Int> best_par_sols(sh.NUM_SEEDS);
    for (int i = 0; i < sh.NUM_SEEDS; ++i){
      best_par_sols[i] = seeds[i];
    }

    std::vector<std::vector<Sol_Int> > sols(sh.NUM_SEEDS);

    qol::CpuTimer record_timer;

    bool flag = false;

    // #pragma omp parallel for shared(flag)// if (sh.NUM_SEEDS > 3)
    for (int seed = 0; (seed < sh.NUM_SEEDS); ++seed){

      if (flag){
        continue;
      }

      // number of merges
      for (int iter = 0; (iter < sh.NUM_ITER) && !flag; ++iter){

        sols[seed].clear();
        sols[seed] = std::vector<Sol_Int>(merge_pop_size);

        double curr_best_obj = best_sol.obj;

        double start_swap = over_timer.elapsedSeconds();

        #pragma omp parallel for shared(curr_best_obj)
        for (int i = 0; i < merge_pop_size; ++i){
          Sol_Int parallel_sol = best_par_sols[seed];
          computeResUse(parallel_sol);

          if (sh.QUIET){
            std::cout <<  "iteration: " << i << " of " << (merge_pop_size * sh.NUM_SEEDS)/omp_get_num_threads() << "\r" << std::flush;
          }
          else{
            std::cout << "full run: " << full_run+1 << " of " << sh.NUM_TOTAL_RUNS << std::endl;
            std::cout << "merge: " << iter << std::endl;
            std::cout << "thread: " << omp_get_thread_num() << std::endl;
            std::cout << "iteration: " << i << " of " << (merge_pop_size * sh.NUM_SEEDS)/omp_get_num_threads() << std::endl << "thread sol: " << seeds[seed].obj <<std::endl;
          }
          // std::cout << "out_streams.size(): " << out_streams.size() << std::endl;
          // (*out_streams[omp_get_thread_num()]) << "printing iteration " << i << " to thread " << omp_get_thread_num() << std::endl;
          // std::cout << "blerp..." << std::endl;
          // ls.SAnoPeriod(parallel_sol);
          // ls.swapWalk(parallel_sol);
          ls.goodSwap(parallel_sol);

          if (sh.OBJ_UB > 0 && parallel_sol.obj > sh.OBJ_UB){
            flag = true;
          }

          computeResUse(parallel_sol);

          curr_best_obj = std::max(curr_best_obj, parallel_sol.obj);

          // std::cout << "blorp" << std::endl;

          if (sh.RECORD_DATA){
            outfile << parallel_sol.obj << std::endl;
            seed_out[seed] << parallel_sol.obj << std::endl;
            seed_out[seed+sh.NUM_SEEDS+1] << seeds[seed].obj << std::endl;
          }

          if (parallel_sol.obj > best_par_sols[seed].obj){
          // test sol
          try{
            bool test_error = verify((*probModel), parallel_sol);
            // bool test_error = false;
            if (!test_error){
              best_par_sols[seed] = parallel_sol;
            }
          } // end try statement
          catch (qol::Exception & ex) {
          }
          catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}
          }

          if (!sh.QUIET){
            std::cout << "objective: " << parallel_sol.obj << std::endl;
            std::cout << "\nbest objective so far (seed " << seed << "): " << best_par_sols[seed].obj << std::endl << std::endl;
          }


          // try{
          //   // std::cout << "\nChecking solution from solver...\n" << std::endl;

          //   bool test_error = verify((*probModel), parallel_sol);

          //   // std::cout << "\nSolution found by merge solver was ";
          //   if (!test_error){
          //     // std::cout << "feasible!\n" << std::endl;
          //   }
          //   else{
          //     // std::cout << "\033[31;1minfeasible!\033[0m\n" << std::endl;
          //     // std::cin.get();
          //     std::cin.get();
          //     throw qol::Exception("Infeasible solution!");
          //     // sol = backup_sol;
          //   }
          // } // end try statement
          // catch (qol::Exception & ex) {
          //   std::cerr << "Error: " << ex.what() << std::endl;
          //   computeResUse(parallel_sol);
          //   ls.repairSolution(parallel_sol);
          // }
          // catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

          // add branch solution to population of solutions
          sols[seed][i] = parallel_sol;
        }

        // output total time taken to generate population
        std::cout << "\ntime taken to generate population: " << over_timer.elapsedSeconds() - start_swap << std::endl << std::endl;

        std::cout << "\n***** Best solution found *****\nseed " << seed << ": " << best_par_sols[seed].obj << std::endl;

        // skip last merge to allow for massive merge
        if (iter < sh.NUM_ITER && !flag){

          std::cout << "************ MERGE " << iter+1 << " OF " << sh.NUM_ITER << " - seed: " << seed << " ************\n" << std::endl;
          qol::CpuTimer merge_timer;
          
          Sol_Int merged_sol;

          doMerge(sm, sols[seed], include, best_par_sols[seed],merged_sol,red_data);

          computeResUse(merged_sol);

          // test final solution for feasibility
          try{
            // std::cout << "\nChecking solution from solver...\n" << std::endl;

            bool test_error = verify((*probModel), merged_sol);

            std::cout << "\nSolution found by merge solver was ";
            if (!test_error){
              std::cout << "feasible!\n" << std::endl;
            }
            else{
              std::cout << "\033[31;1minfeasible!\033[0m\n" << std::endl;
              // std::cin.get();
              throw qol::Exception("Infeasible solution!");
              // sol = backup_sol;
            }
          } // end try statement
          catch (qol::Exception & ex) {
            std::cerr << "Error: " << ex.what() << std::endl;
            computeResUse(merged_sol);
            ls.repairSolution(merged_sol);
          }
          catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

          double iterObj = merged_sol.obj;

          if (sh.OBJ_UB > 0 && merged_sol.obj > sh.OBJ_UB){
            flag = true;
          }

          std::cout << "\nMerge time: " << merge_timer.elapsedSeconds() << std::endl << std::endl;

          Sol_Int sol_iter;
          sol_iter.nT = t_max;

          //TODO: SOMEWHERE HERE IS A SEG FAULT!!! I think it is coming from the fact
          //      that nT for the solutions is not t_max


          double before_timer = over_timer.elapsedSeconds();

          std::cout << "random search value: " << sh.RANDOM_SEARCH << std::endl;
          double after_timer = over_timer.elapsedSeconds();
          // it_times.push_back(after_timer-before_timer);

          std::cout << std::endl << "old best: " << best_sol.obj
                    << std::endl << "branch best: " << best_par_sols[seed].obj
                    << std::endl << "current overall best: " << curr_best_obj
                    << std::endl << "new objective after merge: " << merged_sol.obj;
          if (iterObj > best_par_sols[seed].obj){
              std::cout << std::endl << "***** NEW BEST OBJECTIVE FOUND AFTER MERGE! *****\n" << std::endl;
              bestObj = iterObj;
              best_par_sols[seed] = sol_iter;
          }
          else if (iterObj == bestObj){
              std::cout << std::endl << "***** SAME OBJECTIVE FOUND AFTER MERGE.. *****\n" << std::endl;
          }
          else{
              std::cout << std::endl << "***** WORSE OBJECTIVE FOUND AFTER MERGE! *****\n" << std::endl;
          }
          if (merged_sol.obj > best_par_sols[seed].obj){
            best_par_sols[seed] = merged_sol;
            seeds[seed] = best_par_sols[seed];
          }
        }
      }
    }

    if (flag && sh.RECORD_DATA){
      run_time << run_timer.elapsedSeconds() << std::endl;
      std::cout << "# CPU time        "<<over_timer.elapsedSeconds()
        <<" sec.  Wall:"<< over_timer.elapsedWallTime() <<std::endl;
      continue;
    }

    Sol_Int tot_best_sol = best_par_sols[0];

    std::cout << "Best solutions found on all branches" << std::endl;
    for (int i = 0; i < sh.NUM_SEEDS; ++i){
      if (best_par_sols[i].obj > tot_best_sol.obj){
        tot_best_sol = best_par_sols[i];
      }
      std::cout << "Branch " << i <<  ": " << best_par_sols[i].obj << std::endl;
    }

    // merge best solutions all together
    std::cout << "\nPerforming full merge of all seeds\n" << std::endl;

    for (int p = 0; p < sols.size(); ++p){
      std::cout << "sols[" << p << "]: " << sols[p].size() << std::endl;     
    }
    
    std::vector<Sol_Int> all_sols;
    all_sols.reserve(sh.NUM_SEEDS * sols[0].size());
    for (int seed = 0; seed < sh.NUM_SEEDS; ++seed){
      all_sols.insert(all_sols.end(), sols[seed].begin(), sols[seed].end());
    }




    Sol_Int merged_sol;
    doMerge(sm, all_sols, include, tot_best_sol, merged_sol,red_data);

    // test final solution for feasibility
    try{
      // std::cout << "\nChecking solution from solver...\n" << std::endl;

      bool test_error = verify((*probModel), merged_sol);

      std::cout << "\nSolution found by merge solver was ";
      if (!test_error){
        std::cout << "feasible!\n" << std::endl;
      }
      else{
        std::cout << "\033[31;1minfeasible!\033[0m\n" << std::endl;
        std::cin.get();
        throw qol::Exception("Infeasible solution!");
        // sol = backup_sol;
      }
    } // end try statement
    catch (qol::Exception & ex) {
      std::cerr << "Error: " << ex.what() << std::endl;
      computeResUse(merged_sol);
      ls.repairSolution(merged_sol);
    }
    catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

    std::cout << std::endl << "old best: " << tot_best_sol.obj
              << std::endl << "new objective after merge: " << merged_sol.obj;
    if (merged_sol.obj > tot_best_sol.obj){
        std::cout << std::endl << "***** NEW BEST OBJECTIVE FOUND AFTER MERGE! *****\n" << std::endl;
        tot_best_sol = merged_sol;
    }
    else if (merged_sol.obj == tot_best_sol.obj){
        std::cout << std::endl << "***** SAME OBJECTIVE FOUND AFTER MERGE.. *****\n" << std::endl;
    }
    else{
        std::cout << std::endl << "***** WORSE OBJECTIVE FOUND AFTER MERGE! *****\n" << std::endl;
    }

    double time_now = record_timer.elapsedSeconds();

    if (sh.RECORD_DATA){
      outfile << tot_best_sol.obj << std::endl;
      for (int i = 0; i < sh.NUM_SEEDS; ++i){
        seed_out[i] << tot_best_sol.obj << std::endl;
        seed_out[i+sh.NUM_SEEDS+1] << tot_best_sol.obj << std::endl;
        time_out[i] << time_now << "," << tot_best_sol.obj << std::endl;
      }
    }

    std::cout << "\nPerforming final SA run\n" << std::endl;

    Sol_Int final_sol = tot_best_sol;

    std::vector<Sol_Int> final_sols(merge_pop_size);

    double curr_best_obj = tot_best_sol.obj;

    #pragma omp parallel for
    for (int i = 0; i < merge_pop_size; ++i){
    if (flag){
      continue;
    }
    // for (int i = 0; i < merge_pop_size-1; ++i){
    // for (int i = 0; i < sh.MERGE_POP_SIZE; ++i){
      Sol_Int parallel_sol = final_sol;
      if (sh.QUIET){
        std::cout << "iteration: " << i << " of " << merge_pop_size << "\r" << std::flush;
      }
      else{
        std::cout << "full run: " << full_run+1 << " of " << sh.NUM_TOTAL_RUNS << std::endl;
        std::cout << "thread: " << omp_get_thread_num() << std::endl;
        std::cout << "iteration: " << i << " of " << merge_pop_size << std::endl << "initial value: " << final_sol.obj <<std::endl;
      }
      computeResUse(parallel_sol);
      // ls.SAnoPeriod(parallel_sol);
      ls.swapWalk(parallel_sol);

      if (sh.OBJ_UB > 0 && parallel_sol.obj > sh.OBJ_UB){
        flag = true;
      }

      if (sh.RECORD_DATA){
        outfile << parallel_sol.obj << std::endl;
        seed_out[sh.NUM_SEEDS] << parallel_sol.obj << std::endl;
        seed_out[sh.NUM_SEEDS+sh.NUM_SEEDS+1] << final_sol.obj << std::endl;
        time_out[sh.NUM_SEEDS] << record_timer.elapsedSeconds() << "," << parallel_sol.obj << std::endl;
      }

      if (parallel_sol.obj > tot_best_sol.obj){
        tot_best_sol = parallel_sol;
      }

      if (!sh.QUIET){
        std::cout << "\nbest objective so far:" << tot_best_sol.obj << std::endl << std::endl;
      }

      // add parallel sol to final population
      final_sols[i] = parallel_sol;
    }

    std::cout << "\nPerforming final merge\n" << std::endl;

    if (flag && sh.RECORD_DATA){
      run_time << run_timer.elapsedSeconds() << std::endl;
      std::cout << "# CPU time        "<<over_timer.elapsedSeconds()
        <<" sec.  Wall:"<< over_timer.elapsedWallTime() <<std::endl;
      continue;
    }

    Sol_Int final_merge;
    doMerge(sm, final_sols, include, tot_best_sol, final_merge,red_data);

    // test final solution for feasibility
    try{
      // std::cout << "\nChecking solution from solver...\n" << std::endl;

      bool test_error = verify((*probModel), final_merge);

      std::cout << "\nSolution found by merge solver was ";
      if (!test_error){
        std::cout << "feasible!\n" << std::endl;
      }
      else{
        std::cout << "\033[31;1minfeasible!\033[0m\n" << std::endl;
        // std::cin.get();
        throw qol::Exception("Infeasible solution!");
        // sol = backup_sol;
      }
    } // end try statement
    catch (qol::Exception & ex) {
      std::cerr << "Error: " << ex.what() << std::endl;
      computeResUse(final_merge);
      ls.repairSolution(final_merge);
    }
    catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

    std::cout << std::endl << "old best: " << tot_best_sol.obj
              << std::endl << "new objective after merge: " << final_merge.obj;
    if (final_merge.obj > tot_best_sol.obj){
        std::cout << std::endl << "***** NEW BEST OBJECTIVE FOUND AFTER MERGE! *****\n" << std::endl;
        tot_best_sol = final_merge;
    }
    else if (final_merge.obj == tot_best_sol.obj){
        std::cout << std::endl << "***** SAME OBJECTIVE FOUND AFTER MERGE.. *****\n" << std::endl;
    }
    else{
        std::cout << std::endl << "***** WORSE OBJECTIVE FOUND AFTER MERGE! *****\n" << std::endl;
    }

    if (sh.RECORD_DATA){
      outfile << tot_best_sol.obj << std::endl;
      seed_out[sh.NUM_SEEDS] << tot_best_sol.obj << std::endl;
      seed_out[sh.NUM_SEEDS+sh.NUM_SEEDS+1] << tot_best_sol.obj << std::endl;
      time_out[sh.NUM_SEEDS] << record_timer.elapsedSeconds() << "," << tot_best_sol.obj << std::endl;
      run_time << run_timer.elapsedSeconds() << std::endl;
    }

    std::cout << "Best objective found: " << tot_best_sol.obj << std::endl << std::endl;

      std::cout << "# CPU time        "<<over_timer.elapsedSeconds()
        <<" sec.  Wall:"<< over_timer.elapsedWallTime() <<std::endl;


  // test final solution for feasibility
  try{
    // std::cout << "\nChecking solution from solver...\n" << std::endl;

    bool test_error = verify((*probModel), tot_best_sol);

    std::cout << "\nSolution found by solver was ";
    if (!test_error){
      std::cout << "feasible!\n" << std::endl;
    }
    else{
      std::cout << "\033[31;1minfeasible!\033[0m\n" << std::endl;
      throw qol::Exception("Infeasible solution!");
      // std::cin.get();
      // sol = backup_sol;
    }
  } // end try statement
  catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
  catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

} // end full run

  if (sh.RECORD_DATA){
    outfile.close();
    red_data.close();
    for (int i =0; i < sh.NUM_SEEDS+1; ++i){
      seed_out[i].close();
      seed_out[i+sh.NUM_SEEDS+1].close();
      time_out[i].close();
      run_time.close();
    }
  }


  // free(probModel);
  // return 0;

}

void SinglePSolver::doMerge(SolutionMerger &sm, const std::vector<Sol_Int>&sols, const std::vector<int> &include, Sol_Int &init_sol, Sol_Int &merged_sol, std::ofstream &red_data){
  const int nB = probModel->getNBlock();
  // sm.mergePCPSP(sols,merged);
  std::vector<std::vector<int> > groups;
  std::vector<int> group_map;
  std::vector<int> fixed;
  // sm.mergeCPIT(sols,fixed,groups,group_map);
  // sm.simpleMerge(sols,include,fixed);
  if (sh.MERGE_THRESH){
    sm.fullMergeThresh(sols,include,fixed, groups, group_map);
  }
  else{
    sm.fullMerge(sols,include,fixed, groups, group_map);
  }
  // sols.clear();

  if (sh.RECORD_DATA){
    std::ofstream group_out;
    group_out.open("./tracking_data/" + probModel->getName() + "/GROUPS.csv",std::ios_base::app);
    group_out << groups.size() << std::endl;
    group_out.close();
  }

  // MergeSolverSimple ms = MergeSolverSimple(sh,probModel,init_sol,groups,group_map,fixed);
  MergeSolverCompact ms = MergeSolverCompact(sh,probModel,init_sol,groups,group_map,fixed, include, &red_data);

  ms.solve(merged_sol);
}

// set up MIP and solve
int SinglePSolver::serialMergeSolve(){

  qol::CpuTimer over_timer;

  const int nB = probModel->getNBlock();
  const int t_max = probModel->getNPeriod();
  const int d_max = probModel->getnDestination();
  const int r_max = probModel->getnResources();
  std::vector<Block> * blocks=probModel->getBlock();

  std::ofstream outfile;
  std::ofstream red_data;
  std::vector<std::ofstream> time_out(2);
  std::vector<std::ofstream> restarts_out(1);

  std::ofstream merge_obj;
  merge_obj.open("./sol_files/" + probModel->getName() + "/merge_obj.csv");

  std::vector<std::ofstream> seed_out(4);

  std::ofstream run_time;

  // std::ofstream sol_tracker;
  // sol_tracker.open("sol_tracker.csv");
  // sol_tracker << "SOLUTION TRACKER VALUES" << std::endl;
  // sol_tracker.close();
  //
  // std::ofstream sol_tracker_base;
  // sol_tracker_base.open("sol_tracker_base.csv");
  // sol_tracker_base << "SOLUTION TRACKER BASE VALUES" << std::endl;
  // sol_tracker_base.close();

  if (sh.RECORD_DATA){
  for (int i = 0; i < 1; ++i){
    time_out[i].open("./tracking_data/" + probModel->getName() + "/TIME_"+std::to_string(i)+".csv");
    time_out[i] << "# TIME DATA - SEED " << i << std::endl << "# CPU TIME" << std::endl;
    time_out[i] << "! TIMES" << std::endl;
    time_out[i] << "! LINE 2" << std::endl;

    time_out[i+1].open("./tracking_data/" + probModel->getName() + "/WALL_TIME_"+std::to_string(i)+".csv");
    time_out[i+1] << "# TIME DATA - SEED " << i << std::endl << "# CPU TIME" << std::endl;
    time_out[i+1] << "! WALL_TIMES" << std::endl;
    time_out[i+1] << "! LINE 2" << std::endl;

    restarts_out[i].open("./tracking_data/" + probModel->getName() + "/RESTARTS_"+std::to_string(i)+".csv");
    restarts_out[i] << "# TIME DATA - SEED " << i << std::endl << "# CPU TIME" << std::endl;
    restarts_out[i] << "! RESTART" << std::endl;
    restarts_out[i] << "! LINE 2" << std::endl;
    restarts_out[i] << "0" << std::endl;

    // output files for recording all solutions
    seed_out[i].open("./tracking_data/" + probModel->getName() + "/SEED_"+std::to_string(i)+".csv");
    seed_out[i] << "# PARALLEL MERGE - SEED " << i << std::endl;
    seed_out[i] << "# SA_ALPHA,SA_T_MIN,SA_ITER,FULL_RUNS,MERGE_POP_SIZE,NUM_ITER,ITER_INCR" << std::endl;
    seed_out[i] << "# " << sh.SA_ALPHA <<","<< sh.SA_T_MIN <<","<< sh.SA_ITER <<","<< sh.FULL_RUNS <<","<< sh.MERGE_POP_SIZE <<","<< sh.NUM_ITER << "," << sh.ITER_INCR << std::endl;
    seed_out[i] << "! SOL_IDX " << sh.SOL_IDX << std::endl;
    seed_out[i] << "! LINE 1" << std::endl;
    seed_out[i] << "! NAME " << probModel->getName() << std::endl;
    seed_out[i] << "! POP_SIZE " << sh.MERGE_POP_SIZE << std::endl;
    seed_out[i] << "! SWAPS " << sh.SA_ITER << std::endl;
    seed_out[i] << "! NUM_MERGES " << sh.NUM_ITER << std::endl;
    seed_out[i] << "! MAX_SUB_SWAPS " << sh.MAX_SUB_SWAPS << std::endl;
    seed_out[i] << "! MIP_TIMEOUT " << sh.WINDOW_SEARCH_TIME << std::endl;
    seed_out[i] << "! MIP_REL_GAP " << sh.MIP_GAP << std::endl;


    // output files for recording best solutions
    seed_out[i+1].open("./tracking_data/" + probModel->getName() + "/SEED_"+std::to_string(i)+"_BEST.csv");
    seed_out[i+1] << "# PARALLEL MERGE (BEST SOL)- SEED " << i << std::endl;
    seed_out[i+1] << "# SA_ALPHA,SA_T_MIN,SA_ITER,FULL_RUNS,MERGE_POP_SIZE,NUM_ITER,ITER_INCR" << std::endl;
    seed_out[i+1] << "# " << sh.SA_ALPHA <<","<< sh.SA_T_MIN <<","<< sh.SA_ITER <<","<< sh.FULL_RUNS <<","<< sh.MERGE_POP_SIZE <<","<< sh.NUM_ITER << "," << sh.ITER_INCR << std::endl;
    seed_out[i+1] << "! LINE 2" << std::endl;
    seed_out[i+1] << "! NAME " << probModel->getName() << std::endl;
    seed_out[i+1] << "! POP_SIZE " << sh.MERGE_POP_SIZE << std::endl;
    seed_out[i+1] << "! SWAPS " << sh.SA_ITER << std::endl;
    seed_out[i+1] << "! NUM_MERGES " << sh.NUM_ITER << std::endl;
  }
  run_time.open ("./tracking_data/" + probModel->getName() + "/run_time.csv");
  run_time << "# RUN_TIMES:" << std::endl;
  outfile.open ("./tracking_data/" + probModel->getName() + "/merge_obj.csv");
  red_data.open ("./tracking_data/" + probModel->getName() + "/reduce_data.csv");
  red_data << "# VARIABLES,CONSTRAINTS" << std::endl;
  outfile << "# PARALLEL MERGE" << std::endl;
  outfile << "# SA_ALPHA,SA_T_MIN,SA_ITER,FULL_RUNS,MERGE_POP_SIZE,NUM_ITER,ITER_INCR" << std::endl;
  outfile << "# " << sh.SA_ALPHA <<","<< sh.SA_T_MIN <<","<< sh.SA_ITER <<","<< sh.FULL_RUNS <<","<< sh.MERGE_POP_SIZE <<","<< sh.NUM_ITER << "," << sh.ITER_INCR << std::endl;

  // output files for recording number of groups
  std::ofstream group_out;
  group_out.open("./tracking_data/" + probModel->getName() + "/GROUPS.csv");
  group_out << "# PARALLEL MERGE - " << std::endl;
  group_out << "# SA_ALPHA,SA_T_MIN,SA_ITER,FULL_RUNS,MERGE_POP_SIZE,NUM_ITER,ITER_INCR" << std::endl;
  group_out << "# " << sh.SA_ALPHA <<","<< sh.SA_T_MIN <<","<< sh.SA_ITER <<","<< sh.FULL_RUNS <<","<< sh.MERGE_POP_SIZE <<","<< sh.NUM_ITER << "," << sh.ITER_INCR << std::endl;
  group_out << "! LINE 2" << std::endl;
  group_out << "! NAME " << probModel->getName() << std::endl;
  group_out << "! POP_SIZE " << sh.MERGE_POP_SIZE << std::endl;
  group_out << "! GROUPS " << std::endl;
  group_out.close();
  }

  BranchNode_info base_info(nB,t_max,d_max,r_max);
  std::vector<int> include(nB,0);

  computeUPIT(base_info, include);

  for (int full_run = 0; full_run < sh.NUM_TOTAL_RUNS; ++full_run){

    qol::CpuTimer run_timer;

    std::vector<Sol_Int> seeds;

    int merge_pop_size = sh.MERGE_POP_SIZE;

    int n_periods = t_max;

    std::string lpFile="";// ,solnFile = probModel->getName()+".sol";

    SolutionMerger sm(sh);
    LocalSearch ls = LocalSearch(sh,probModel,include);

    std::cout << probModel->getName() << std::endl;

    if (sh.LOAD_SEEDS){
      loadSols(seeds);
    }
    else{
      std::cout << "getting seeds.. "<<std::flush;
      getSeeds(base_info, seeds);
      std::cout << "done!" << std::endl;

      // std::string path_name = "./rand_sols/";

      // saveSols(seeds, path_name);
      // return;
    }

    Sol_Int best_sol = seeds[0];

    // Sol_Int merged_sol;

    // doMerge(sm, seeds, include, seeds[0] ,merged_sol,red_data);

    double bestObj = best_sol.obj;


    std::vector<Sol_Int> sols;
    qol::CpuTimer record_timer;

    bool flag = false;

    int resets = 0;

    // number of merges
    for (int iter = 0; (iter < sh.NUM_ITER) && !flag; ++iter){
      std::vector<std::vector<double> > data_out(4, std::vector<double>(merge_pop_size));
      std::vector<int> idx_order(merge_pop_size);
      int order_count = 0;

      sols.clear();
      sols = std::vector<Sol_Int>(merge_pop_size);

      double curr_best_obj = best_sol.obj;

      double start_swap = over_timer.elapsedSeconds();

      std::cout << "seed solution: " << seeds[0].obj << std::endl;

      // save seed solution
      sols[0] = seeds[0];

      data_out[0][0] = sols[0].obj;
      data_out[1][0] = best_sol.obj;
      data_out[2][0] = record_timer.elapsedSeconds();
      data_out[3][0] = record_timer.elapsedWallTime();

      // write out current solution
      std::string path_name = "./sol_files/";
      std::vector<Sol_Int> temp_vec;
      temp_vec.push_back(seeds[0]);
      saveSols(temp_vec, path_name, iter);
      merge_obj << iter << " " << seeds[0].obj << std::endl;

      #pragma omp parallel for shared(best_sol,order_count)
      for (int i = 1; i < merge_pop_size; ++i){
        // Sol_Int sol = best_sol;
        Sol_Int sol = seeds[0];
        computeResUse(sol);

        if (sh.QUIET){
          if (omp_get_thread_num() == 0){
            std::cout <<  "iteration: " << i*omp_get_num_threads() << " of " << merge_pop_size << "\r" << std::flush;
          }
        }
        else{
          std::cout << "full run: " << full_run+1 << " of " << sh.NUM_TOTAL_RUNS << std::endl;
          std::cout << "merge: " << iter << std::endl;
          std::cout << "thread: " << omp_get_thread_num() << std::endl;
          std::cout << "iteration: " << i << " of " << (merge_pop_size)/omp_get_num_threads() << std::endl << "thread sol: " << seeds[0].obj <<std::endl;
        }

        // ls.swapWalk(sol);
        ls.goodSwap(sol);

        if (sh.OBJ_UB > 0 && sol.obj > sh.OBJ_UB){
          flag = true;
        }

        computeResUse(sol);

        if (sol.obj > best_sol.obj){
          // test sol
          try{
            bool test_error = verify((*probModel), sol);
            // bool test_error = false;
            if (!test_error){
              best_sol = sol;
            }
          } // end try statement
          catch (qol::Exception & ex) {
          }
          catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}
        }

        if (sh.RECORD_DATA){
          data_out[0][i] = sol.obj;
          data_out[1][i] = best_sol.obj;
          data_out[2][i] = record_timer.elapsedSeconds();
          data_out[3][i] = record_timer.elapsedWallTime();
          idx_order[i] = order_count++;
        }

        if (!sh.QUIET){
          std::cout << "objective: " << sol.obj << std::endl;
          std::cout << "\nbest objective so far: " << best_sol.obj << std::endl << std::endl;
        }


        // try{
        //   // std::cout << "\nChecking solution from solver...\n" << std::endl;
        //
        //   bool test_error = verify((*probModel), sol);
        //
        //   // std::cout << "\nSolution found by merge solver was ";
        //   if (!test_error){
        //     // std::cout << "feasible!\n" << std::endl;
        //   }
        //   else{
        //     // std::cout << "\033[31;1minfeasible!\033[0m\n" << std::endl;
        //     // std::cin.get();
        //     std::cin.get();
        //     throw qol::Exception("Infeasible solution!");
        //     // sol = backup_sol;
        //   }
        // } // end try statement
        // catch (qol::Exception & ex) {
        //   std::cerr << "Error: " << ex.what() << std::endl;
        //   computeResUse(sol);
        //   ls.repairSolution(sol);
        // }
        // catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}


        // add branch solution to population of solutions
        sols[i] = sol;
      }
      
      if (flag){
        break;
      }

      if (sh.RECORD_DATA){
        // sort data by order_count
        std::priority_queue<std::pair<int, int> > q;
        for (int i = 0; i < data_out[0].size(); ++i) {
          // q.push(std::pair<double, int >(-data_out[2][i], i));
          q.push(std::pair<double, int >(-idx_order[i], i));
        }

        double best_obj = -1e20;

        while(!q.empty()) {
          int data_idx = q.top().second;
          seed_out[0] << data_out[0][data_idx] << std::endl;
          if (data_out[1][data_idx] > best_obj){
            seed_out[1] << data_out[1][data_idx] << std::endl;
            best_obj = data_out[1][data_idx];
          }
          else{
            seed_out[1] << best_obj << std::endl;
          }
          time_out[0] << data_out[2][data_idx] << std::endl;
          time_out[1] << data_out[3][data_idx] << std::endl;
          q.pop();
        }
      }

      // output total time taken to generate population
      std::cout << "\ntime taken to generate population: " << over_timer.elapsedSeconds() - start_swap << std::endl << std::endl;

      std::cout << "\n***** Best solution found *****\n" << best_sol.obj << std::endl;

      std::cout << "************ " << probModel->getName() << " MERGE " << iter+1 << " OF " << sh.NUM_ITER << " (resets: " << resets << ") ************\n" << std::endl;
      qol::CpuTimer merge_timer;

      Sol_Int merged_sol(nB,t_max);

      // TODO: MAYBE INIT MERGED SOL!

      doMerge(sm, sols, include, best_sol, merged_sol, red_data);

      computeResUse(merged_sol);

      // test final solution for feasibility
      try{
        // std::cout << "\nChecking solution from solver...\n" << std::endl;

        bool test_error = verify((*probModel), merged_sol);

        std::cout << "\nSolution found by merge solver was ";
        if (!test_error){
          std::cout << "feasible!\n" << std::endl;
        }
        else{
          std::cout << "\033[31;1minfeasible!\033[0m\n" << std::endl;
          // std::cin.get();
          throw qol::Exception("Infeasible solution!");
          // sol = backup_sol;
        }
      } // end try statement
      catch (qol::Exception & ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        computeResUse(merged_sol);
        ls.repairSolution(merged_sol);

        // iter--;
        resets++;
        restarts_out[0] << resets << std::endl;
      }
      catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

      double iterObj = merged_sol.obj;

      if (sh.OBJ_UB > 0 && merged_sol.obj > sh.OBJ_UB){
        flag = true;
      }

      std::cout << "\nMerge time: " << merge_timer.elapsedSeconds() << std::endl << std::endl;

      Sol_Int sol_iter;
      sol_iter.nT = t_max;

      //TODO: SOMEWHERE HERE IS A SEG FAULT!!! I think it is coming from the fact
      //      that nT for the solutions is not t_max


      double before_timer = over_timer.elapsedSeconds();

      std::cout << "random search value: " << sh.RANDOM_SEARCH << std::endl;
      double after_timer = over_timer.elapsedSeconds();
      // it_times.push_back(after_timer-before_timer);

      std::cout << std::endl << "old best: " << best_sol.obj
                << std::endl << "new objective after merge: " << merged_sol.obj;
      if (merged_sol.obj > best_sol.obj){
          std::cout << std::endl << "***** NEW BEST OBJECTIVE FOUND AFTER MERGE! *****\n" << std::endl;
          best_sol = merged_sol;
          seeds[0] = best_sol;
      }
      else if (merged_sol.obj == best_sol.obj){
          std::cout << std::endl << "***** SAME OBJECTIVE FOUND AFTER MERGE.. *****\n" << std::endl;
      }
      else{
          std::cout << std::endl << "***** WORSE OBJECTIVE FOUND AFTER MERGE! *****\n" << std::endl;

          try{
            // std::cout << "\nChecking solution from solver...\n" << std::endl;

            bool test_error = verify((*probModel), best_sol);

            // std::cout << "\nSolution found by merge solver was ";
            if (!test_error){
              // std::cout << "feasible!\n" << std::endl;
              seeds[0] = best_sol;
            }
          } // end try statement
          catch (qol::Exception & ex) {
            // std::cerr << "Error: " << ex.what() << std::endl;
            // computeResUse(merged_sol);
            // ls.repairSolution(merged_sol);

            // iter--;
            // resets++;
            // restarts_out[0] << resets << std::endl;
          }
          catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}
      }
    }

    if (sh.RECORD_DATA){
      outfile << best_sol.obj << std::endl;
      for (int i = 0; i < 2; ++i){
      seed_out[i] << best_sol.obj << std::endl << "! RESTARTS " << resets << std::endl
                  << "! CPU_TIME " << over_timer.elapsedSeconds() << std::endl
                  << "! WALL_TIME " << over_timer.elapsedWallTime() << std::endl;
      }
      time_out[0] << record_timer.elapsedSeconds() << std::endl;
      time_out[1] << record_timer.elapsedWallTime() << std::endl;
      run_time << run_timer.elapsedSeconds() << std::endl;
    }

    std::cout << "Best objective found: " << best_sol.obj << std::endl << std::endl;

      std::cout << "# CPU time        "<<over_timer.elapsedSeconds()
        <<" sec.  Wall:"<< over_timer.elapsedWallTime() <<std::endl;

      std::cout <<"\nrun time:" << std::endl;
      std::cout << "# CPU time        "<<run_timer.elapsedSeconds()
        <<" sec.  Wall:"<< run_timer.elapsedWallTime() <<std::endl;


} // end full run

  if (sh.RECORD_DATA){
    outfile.close();
    red_data.close();
    for (int i =0; i < 2; ++i){
      seed_out[i].close();
      seed_out[i+2].close();
      time_out[i].close();
      restarts_out[i].close();
      run_time.close();
    }
  }


  // free(probModel);
  // return 0;

}

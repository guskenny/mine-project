#include "LocalSearch.h"

const int LocalSearch::getDirection(){
    std::uniform_int_distribution<int> uni(0,NUM_DIRS-1);
    return uni(rng);
}

bool LocalSearch::computeCone(const Sol_Int &sol, const int &cone_tip,
  SetObj &cone, long double &cone_profit, std::vector<long double> &cone_res_use,
  const int &direction, int &cone_depth){

  const int nB = model->graph.getNumNodes();
  const int r_max = model->getnResources();
  const int t_max = model->getNPeriod();

  // reset all data structures
  cone.clear();
  cone_res_use = std::vector<long double>(r_max,0.0);
  cone_profit = 0.0;
  cone_depth = 1;

  // create stack for cone search
  std::vector<std::pair<int,int> > stack;
  stack.push_back(std::make_pair(cone_tip,1));
  cone.addElement(cone_tip);

  while (!stack.empty()){
    int block_idx = stack.back().first; // get next block id
    int depth = stack.back().second;
    cone_depth = std::max(cone_depth, depth);
    const Block & block=model->getBlock(block_idx);
    stack.pop_back(); // remove block from stack

    cone_profit += block.getProfit(0); // function returns un-discounted

    for (int r=0; r < r_max; ++r){
      long double block_res = block.getRCoef(0,r);
      cone_res_use[r] += block_res;
    }

    // get graph node of current block
    const Node* curr = model->graph.getNode(block_idx);

    std::vector<Node*> connected = std::vector<Node*>();

    const int num_connected = curr->getConnected(direction, connected);

    // iterate over all connected nodes to current block
    for (size_t node = 0; node < num_connected; ++node){
      int node_idx = connected[node]->getID();
      // if successor in same period and not already added (or not in include), add to stack
      if (sol.x[node_idx] == sol.x[cone_tip]){
        if (!cone.is_element(node_idx) && include[node_idx]){
          cone.addElement(node_idx);
          stack.push_back(std::make_pair(node_idx, depth+1));
        }
      }
    }
  }

  return true;
}

void LocalSearch::swapConePeriod(Sol_Int &sol, const int period, const int dir,
  const SetObj &cone, const long double cone_profit, const std::vector<long double> &cone_res_use, std::vector<SetObj> &period_blocks){

    const int nB = model->graph.getNumNodes();
    const int r_max = model->getnResources();
    const double rate = model->getDiscountRate();

    int swap_period = period + dir - !dir;

    // for every block in the cone, swap to swap period
    for (auto it = cone.set_data.begin(); it != cone.set_data.end(); ++it){
      // if (!period_blocks[period].is_element(*it)){
      //   std::cout << "error: block " << *it << " not in period_blocks["<< period<< "]!!" << std::endl;
      //   std::cin.get();
      // }
      // if (period_blocks[swap_period].is_element(*it)){
      //   std::cout << "error: block " << *it << " already in period_blocks["<< swap_period<< "]!!" << std::endl;
      //   std::cin.get();
      // }

      // std::cout << "swapping block " << *it << " from " << period << " to " << swap_period << std::endl;
      period_blocks[period].removeElement(*it);
      period_blocks[swap_period].addElement(*it);
      sol.x[*it] = swap_period;
    }

    // std::cout << "testing period_blocks... " << std::flush;
    // for (int b = 0; b < nB; ++b){
    //   if (!period_blocks[sol.x[b]].is_element(b)){
    //     std::cout << "error: block " << b << " not in period_blocks["<< sol.x[b] << "]!!" << std::endl;
    //     std::cin.get();
    //   }
    // }
    // std::cout << "done!" << std::endl;

    // update profit
    sol.obj += (cone_profit / pow(1+rate,swap_period)) -
      (cone_profit / pow(1+rate,period));

    // update resource usage
    for (int r = 0; r < r_max; ++r){
      sol.res_use[period][r] -= cone_res_use[r];
      sol.res_use[swap_period][r] += cone_res_use[r];
    }

}

// TODO: FIX REPAIR SOLUTION, NOT DOING ANY SWAPS! MAYBE UPDATE RES_USE

void LocalSearch::repairSolution(Sol_Int &sol){
  const int nB = model->graph.getNumNodes();
  const int r_max = model->getnResources();
  const int t_max = model->getNPeriod();

  std::cout << "repairing solution... " << std::flush;

  bool flag = true;
  int fix_count = 0;
  while (flag && fix_count < 100){
    fix_count++;
    flag = false;
    // verify precedence constraints
  	for(int a=0; a<nB; a++){
  		if(sol.x[a]<t_max){
  			const Block & block_a = model->getBlock(a); //successor
  			const std::vector<int> & preds = block_a.getPreds();
  			const int n = preds.size();
  			for(int p=0; p<n; p++){
  				const int b = preds[p];    // predecessor
  				if(sol.x[a] < sol.x[b] && sol.x[b] != t_max){
            for (int r = 0; r < r_max; ++r){
              long double block_res = block_a.getRCoef(0,r);
              sol.res_use[sol.x[a]][r] -= block_res;
              sol.res_use[sol.x[b]][r] += block_res;
            }
            sol.x[a] = sol.x[b];
            flag = true;
  				}
  			}
  		}
  	}
  }

  std::vector<SetObj> period_blocks(t_max+1, SetObj(nB));
  for (int b = 0; b < nB; ++b){
      period_blocks[sol.x[b]].addElement(b);
  }

  SetObj cone(nB);
  int cone_depth = 0;
  std::vector<long double> cone_res_use(NUM_DIRS,0.0);
  long double cone_profit = 0.0;

  int swap_period = -1;

  std::vector<double> res_limits(r_max,0.0);
  for (size_t r = 0;r < r_max; ++r){
    res_limits[r] = model->getLimit(r, 0);
  }

  Sol_Int backup_sol = sol;
  std::vector<SetObj> backup_period_blocks = period_blocks;
  bool skip = false;

  int skips = 0;

  // iterate over periods and push blocks later until resource limits are fine
  for (int period = 0; period < t_max; ++period){
    skip = true;
    for (int r = 0; r < r_max; ++r){
      if (sol.res_use[period][r] > res_limits[r]){
        skip = false;
      }
    }

    int dir = FORWARD;
    std::vector<int> boundary(nB,0);
    while(!skip){
      // std::cout << "swapping!" << std::endl;
      int cone_tip = period_blocks[period].getRandomElement(rng);

      int best_tip = cone_tip;
      long double best_profit = -1e50;

      for (int s = 0; s < sh.SWAP_POP; ++s){
        // computing cone from cone tip in direction dir
        computeCone(sol, cone_tip, cone, cone_profit, cone_res_use,
          dir, cone_depth);
        if ((cone_profit/cone.get_set_size()) > best_profit){
          best_profit = (cone_profit/cone.get_set_size());
          // std::cout << "updating best profit to "<< best_profit << "!" << std::endl;
          best_tip = cone_tip;
        }
        cone_tip = period_blocks[period].getRandomElement(rng);
        // std::cout << "new cone tip: " << cone_tip << ", period: " << sol.x[cone_tip] << std::endl;
      }

      // computing best cone from cone tip in direction dir
      computeCone(sol, best_tip, cone, cone_profit, cone_res_use,
        dir, cone_depth);

      // swap blocks from current period to swap period and update resource usage
      // profit and period_blocks
      swapConePeriod(sol, period, dir, cone, cone_profit, cone_res_use, period_blocks);

      skip = true;
      for (int r = 0; r < r_max; ++r){
        if (sol.res_use[period][r] > res_limits[r]){
          skip = false;
        }
      }
    }
  }

  std::cout << "done!" << std::endl;

  std::cout << "testing solution new solution... " << std::flush;
  try{
    // std::cout << "\nChecking solution from solver...\n" << std::endl;

    bool test_error = verify((*model), sol);

    // std::cout << "Solution found by solver was ";
    if (!test_error){
      std::cout << "feasible!\n" << std::endl;
      // std::cin.get();
    }
    else{
      std::cout << "\033[31;1minfeasible!\033[0m\n" << std::endl;
      // std::cin.get();
      throw qol::Exception("Infeasible solution!");
    }
  } // end try statement
  catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
  catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

}

void LocalSearch::swapWalk(Sol_Int &sol){
  const int nB = model->graph.getNumNodes();
  const int r_max = model->getnResources();
  const int t_max = model->getNPeriod();

  // std::ofstream sol_tracker;
  // sol_tracker.open("sol_tracker.csv", std::ios_base::app);
  //
  // std::ofstream sol_tracker_base;
  // sol_tracker_base.open("sol_tracker_base.csv", std::ios_base::app);

  double base_obj = sol.obj;

  // set up data structures
  SetObj cone(nB);
  if (!sh.QUIET){
    std::cout << "Beginning swapWalk:" << std::endl;
  }
  std::vector<SetObj> period_blocks(t_max+1, SetObj(nB));
  for (int b = 0; b < nB; ++b){
      period_blocks[sol.x[b]].addElement(b);
  }

  // std::cout << "testing period_blocks... " << std::flush;
  for (int b = 0; b < nB; ++b){
    if (!period_blocks[sol.x[b]].is_element(b)){
      std::cout << "error: block " << b << " not in period_blocks["<< sol.x[b] << "]!!" << std::endl;
      std::cin.get();
    }
  }
  // std::cout << "done!" << std::endl;

  std::vector<long double> cone_res_use(NUM_DIRS,0.0);
  long double cone_profit = 0.0;
  int swap_period = -1;

  std::vector<double> res_limits(r_max,0.0);

  for (size_t r = 0;r < r_max; ++r){
    res_limits[r] = model->getLimit(r, 0);
  }

  double best_obj = 0.0;
  int cone_depth = 0;

  bool skip = false;

  int skips = 0;

  // iterate for number of swaps
  for (int swap = 0; swap < sh.SA_ITER && sol.obj < sh.OBJ_UB; ++swap){
    for (int sub_swap = 0; sub_swap < sh.FULL_RUNS && sol.obj < sh.OBJ_UB; ++sub_swap){
      Sol_Int backup_sol = sol;
      std::vector<SetObj> backup_period_blocks = period_blocks;
      // if (swap % 3 == 0){
      //   sol_tracker << sol.obj << std::endl;
      //   sol_tracker_base << base_obj << std::endl;
      // }
      // get random block for cone tip
      std::uniform_int_distribution<int> uni(0,nB-1);
      int cone_tip = uni(rng);
      // if invalid block, select again
      while (sol.x[cone_tip] == t_max || !include[cone_tip]){
        cone_tip = uni(rng);
      }

      int period = sol.x[cone_tip];
      std::vector<int> boundary(nB,0);
      int dir = getDirection();

      // if cone in period 0, direction forced forwards
      if (period == 0){
        dir = FORWARD;
      }
      if (period == t_max-1){
        dir = BACKWARD;
      }

      int best_tip = cone_tip;
      // best_profit is direction dependent, if backwards we want to move the blocks with the least profit density
      long double best_profit = -1e50;

      for (int s = 0; s < sh.SWAP_POP; ++s){
        // computing cone from cone tip in direction dir
        computeCone(sol, cone_tip, cone, cone_profit, cone_res_use,
          dir, cone_depth);
        if ((cone_profit/cone.get_set_size()) > best_profit){
          best_profit = (cone_profit/cone.get_set_size());
          // std::cout << "updating best profit to "<< best_profit << "!" << std::endl;
          best_tip = cone_tip;
        }
        cone_tip = period_blocks[period].getRandomElement(rng);
        // std::cout << "new cone tip: " << cone_tip << ", period: " << sol.x[cone_tip] << std::endl;
      }

      // computing cone from cone tip in direction dir
      computeCone(sol, best_tip, cone, cone_profit, cone_res_use,
        dir, cone_depth);

      // swap blocks from current period to swap period and update resource usage
      // profit and period_blocks
      swapConePeriod(sol, period, dir, cone, cone_profit, cone_res_use, period_blocks);

      // create vector with both current swap periods in
      std::vector<int> swap_periods = {period, period + dir - !dir};

      if (sh.AUTO_REPAIR){
        repairSolution(sol);
      }
      else{
        int sub_sub_swaps = 1;
        bool flag = false;
        // keep swapping while resources dont match up
        while(!flag){
          flag = true;
          int temp_t = 0;
          // check all resources to see if they are over the limit
          for (int r = 0; r < r_max && flag; ++r){
            for (int t = 0; t < swap_periods.size() && flag; ++t){
              if (sol.res_use[swap_periods[t]][r] > model->getLimit(r, swap_periods[t])){
                flag = false;
                temp_t = t;
                break;
              }
            }
          }
          if (!flag){
            boundary = std::vector<int>(nB,0);
            int temp_dir = dir != temp_t; // if t is 0 then temp_dir = dir, else temp_dir = !dir
            int temp_p = swap_periods[temp_t];
            // select random block from temp_p
            int temp_b = period_blocks[temp_p].getRandomElement(rng);

            // SetObj temp_group(nB);
            // for(int b = 0; b < nB; ++b){
            //   if (sol.x[b] == temp_p)
            //     temp_group.addElement(b);
            // }
            // int temp_b = temp_group.getRandomElement(rng);

            // // testing swapping beause its going haywire
            // std::cout << "random block " << temp_b << " selected from period " << temp_p << ", block actually in " << sol.x[temp_b] << " [" << swap_periods[0] << "," << swap_periods[1] << "]" << std::flush;
            //
            // if (temp_p != sol.x[temp_b]){
            //     std::cout << ", block";
            //     if (!period_blocks[temp_p].is_element(temp_b))
            //       std::cout << " NOT";
            //     std::cout << " in period_blocks[" << temp_p << "], and";
            //     if (!period_blocks[sol.x[temp_b]].is_element(temp_b))
            //       std::cout << " NOT";
            //     std::cout << " in period_blocks[" << sol.x[temp_b] << "] ";
            //     std::cout << "\033[31;1m*** ERROR ***\033[0m\n";
            //     std::cin.get();
            // }
            // std::cout << std::endl;

            best_tip = temp_b;
            best_profit = -1e50;

            for (int s = 0; s < sh.SWAP_POP && best_profit < 0; ++s){
              // computing cone from cone tip in direction dir
              computeCone(sol, temp_b, cone, cone_profit, cone_res_use,
                temp_dir, cone_depth);
              if ((cone_profit/cone.get_set_size()) > best_profit){
                best_profit = (cone_profit/cone.get_set_size());
                // std::cout << "updating best profit to "<< best_profit << "!" << std::endl;
                best_tip = temp_b;
              }
              temp_b = period_blocks[temp_p].getRandomElement(rng);
              // std::cout << "new cone tip: " << cone_tip << ", period: " << sol.x[cone_tip] << std::endl;
            }


            // compute cone in temp_dir direction
            computeCone(sol, best_tip, cone, cone_profit, cone_res_use,
              temp_dir, cone_depth);
            // swap blocks from current period to swap period and update resource usage
            // profit and period_blocks
            swapConePeriod(sol, temp_p, temp_dir, cone, cone_profit, cone_res_use, period_blocks);
            // std::cout << "swapping from " << swap_periods[0] << " to " << swap_periods[1] << ", current objective value: " << sol.obj << std::endl;

            sub_sub_swaps++;
            if (sub_sub_swaps > sh.MAX_SUB_SWAPS){
              skips++;
              // std::cout << "not frozen, just swapping from " << temp_p << " to " << int(temp_p + temp_dir - !temp_dir) << ", current objective value: " << sol.obj << std::endl;
              // std::cout << "sub_sub_swaps over 100, skipping current sub_swap!" << std::endl;
              // drop current sub-swap
              // swap--;
              sol = backup_sol;
              period_blocks = backup_period_blocks;
              flag = true;
              // std::cin.get();
            }
          }
        }
      }
      // std::cout << "swapped from period " << swap_periods[0] << " to " << swap_periods[1] << std::endl;
      // if (sub_sub_swaps > 1)
      //   std::cout << "\033[32;1m";
      // std::cout << "swap: " << swap << ", sub_swap: "<< sub_swap << ", " << sub_sub_swaps << " sub-sub-swaps made, current objective value: " << sol.obj << "\033[0m" << std::endl;

      // // make new backup current sol
      // try{
      //   // std::cout << "\nChecking solution from solver...\n" << std::endl;
      //
      //   bool test_error = verify((*model), sol);
      //
      //   // std::cout << "Solution found by solver was ";
      //   if (!test_error){
      //     // std::cout << "feasible!\n" << std::endl;
      //     // backup_sol = sol;
      //     // backup_period_blocks = period_blocks;
      //   }
      //   else{
      //     std::cout << "\033[31;1minfeasible!\033[0m\n" << std::endl;
      //     std::cin.get();
      //     throw qol::Exception("Infeasible solution!");
      //     // sol = backup_sol;
      //   }
      // } // end try statement
      // catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
      // catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

    }
  }
  if (!sh.QUIET){
    std::cout << "number of skips: " << skips << std::endl;
  }
  // sol_tracker.close();
  // sol_tracker_base.close();
}

// void LocalSearch::SAnoPeriod(Sol_Int &sol){
//   const int nB = model->graph.getNumNodes();
//   const int r_max = model->getnResources();
//   const int t_max = model->getNPeriod();
//
//   // set up data structures
//   std::vector<SetObj> bounds(NUM_DIRS, SetObj(nB));
//   SetObj cone(nB);
//
//   std::vector<double> cone_res_use(NUM_DIRS,0.0);
//   double cone_profit = 0.0;
//   int swap_period = -1;
//
//   std::vector<double> res_limits(r_max,0.0);
//
//   for (size_t r = 0;r < r_max; ++r){
//     res_limits[r] = model->getLimit(r, 0);
//   }
//
//   double best_obj = 0.0;
//   int best_obj_run = 0;
//
//   int abs_max_cone_depth = 0;
//   int abs_max_cone_size = 0;
//
//   std::cout << "Beginning SAnoPeriod:" << std::endl;
//
//   for (int run_no = 0; run_no < sh.FULL_RUNS && sol.obj < sh.OBJ_UB; ++run_no){
//       // temperature for simulated annealing
//       double temp = 1.0;
//       double old_mine_total = sol.obj;
//
//       std::cout << "Run: " << run_no+1 << "/" << sh.FULL_RUNS <<": " << old_mine_total << std::endl;
//
//       bool no_boundary = false;
//       int num_swaps = 0;
//       int num_iters = 0;
//
//       int max_cone_size = 0;
//       int tot_cone_size = 0;
//       int max_cone_depth = 0;
//       int tot_cone_depth = 0;
//       int cone_depth = 0;
//
//       bool period_swap = false;
//
//       while (temp > sh.SA_T_MIN && sol.obj < sh.OBJ_UB){
//         for (int i = 0; i < sh.SA_ITER; ++i){
//
//           // std::cout << "getting random direction... " << std::flush;
//           // get random direction
//           const int dir = getDirection();
//           // std::cout << "done!" << std::endl;
//
//           // std::cout << "getting random block... " << std::flush;
//           // get random block for cone tip
//           std::uniform_int_distribution<int> uni(0,nB-1);
//           int cone_tip = uni(rng);
//           while (sol.x[cone_tip] == t_max){
//             // std::cout << "invalid cone" << std::endl;
//             cone_tip = uni(rng);
//           }
//           // std::cout << "done!" << std::endl;
//
//           // std::cout << "computing cone... " << std::flush;
//           // compute cone
//
//           // if (sh.BOUND_DEPTH == 1){
//           //   computeSingleBlock(sol, cone_tip, cone, cone_profit, cone_res_use, swap_period, dir);
//           // }
//           // else{
//
//             // computeCone(sol, cone_tip, cone, cone_profit, cone_res_use,
//               // swap_period, dir, cone_depth);
//             // if (dir == FORWARD){
//             //   swap_period = sol.x[cone_tip]+1;
//             // }
//             // else{
//             //   if (sol.x[cone_tip] > 0){
//             //     swap_period = sol.x[cone_tip]-1;
//             //   }
//             // }
//
//             // if (swap_period > sol.x[cone_tip]+1 || swap_period < sol.x[cone_tip]-1){
//             //   std::cout << "shifted " << sol.x[cone_tip] << " to " << swap_period << "!" << std::endl;
//             // }
//             // std::cout << swap_period << " " << std::flush;
//           // }
//
//           // std::cout << "done!" << std::endl;
//
//           // std::cout << "deciding move... " << std::flush;
//           // decide if move is going to be taken
//           // std::cout << "test 3..." << std::endl;
//           if (decideMove(sol, swap_period, cone_tip, cone_res_use,
//             res_limits, cone_profit, temp)){
//               // std::cout << "test 3.9...(true)" << std::endl;
//               // std::cout << "accepted!" << std::endl;
//               num_swaps++;
//               // swapConePeriod(sol, swap_period, cone);
//               period_swap = true;
//               max_cone_size = std::max(max_cone_size, cone.get_set_size());
//               tot_cone_size += cone.get_set_size();
//               max_cone_depth = std::max(max_cone_depth, cone_depth);
//               tot_cone_depth += cone_depth;
//           }
//           else{
//             // std::cout << "test 3.9(false)..." << std::endl;
//             period_swap = false;
//             // std::cout << "rejected!" << std::endl;
//           }
//           cone.clear();
//           num_iters++;
//         }
//         // reduce temperature
//         temp *= sh.SA_ALPHA;
//       }
//       std::cout << "Number of swaps made: " << num_swaps
//                 << " / " << num_iters << "\nNew total: ";
//       if (sol.obj > old_mine_total){
//         std::cout << "\e[32m\e[1m";
//         best_obj = sol.obj;
//         best_obj_run = run_no+1;
//       }
//       else if (sol.obj < old_mine_total){
//         std::cout << "\e[31m\e[1m";
//       }
//       std::cout << sol.obj << "\e[0m" << std::endl
//                 << "avg cone size: " << (double)tot_cone_size/num_swaps << ", max cone size: " << max_cone_size << std::endl
//                 << "avg cone depth: " << (double)tot_cone_depth/num_swaps << ", max cone depth: " << max_cone_depth << std::endl;
//       abs_max_cone_depth = std::max(abs_max_cone_depth,max_cone_depth);
//       abs_max_cone_size = std::max(abs_max_cone_size,max_cone_size);
//
//   }
//
//   std::cout << "\nBest solution found: " << best_obj << ", found at run: " << best_obj_run << std::endl << std::endl;
//
//   try{
//     std::cout << "\nChecking solution from solver...\n" << std::endl;
//
//     bool test_error = verify((*model), sol);
//
//     std::cout << "Solution found by solver was ";
//     if (!test_error)
//       std::cout << "feasible!\n" << std::endl;
//     else{
//       std::cout << "infeasible!\n" << std::endl;
//       throw qol::Exception("Infeasible solution!");
//     }
//   } // end try statement
//   catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
//   catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}
//
//   std::cout << "\ntotal max cone size: " << abs_max_cone_size << std::endl;
//   std::cout << "total max cone depth: " << abs_max_cone_depth << std::endl <<std::endl;
//
//   //std::cin.get();
// }

// void LocalSearch::runSingleTests(Sol_Int &sol){
//
//   sol.x = std::vector<int>{1,0,0,0,0,0,0,0,0,0,0,
//                            1,1,0,0,0,0,0,0,0,0,0,
//                            1,1,1,3,0,0,0,0,0,0,2,
//                            1,1,3,3,3,0,0,0,0,2,2,
//                            1,3,3,3,3,3,0,2,2,2,2,
//                            5,3,3,3,3,3,4,2,2,2,2,
//                            5,5,3,3,3,4,4,4,2,2,2,
//                            5,5,5,5,4,4,4,4,4,4,2,
//                            5,5,5,5,5,4,4,4,4,4,4,
//                            5,5,5,5,5,5,4,4,4,4,4,
//                            5,5,5,5,5,5,5,4,4,4,4};
//
//   int NUM_BLOCKS = sol.x.size();
//   int NUM_PERIODS = *max_element(sol.x.begin(), sol.x.end())+1;
//   int WIDTH = (int)sqrt(NUM_BLOCKS);
//   Graph graph = Graph();
//
//   genTestGraph(NUM_BLOCKS,WIDTH,graph);
//
//   model->graph = graph;
//
//   // for (int b = 0; b < NUM_BLOCKS; ++b){
//   //   if (b % WIDTH == 0){
//   //     std::cout << std::endl << std::endl;
//   //   }
//   //
//   //   const Node* curr = graph.getNode(b);
//   //
//   //   std::cout << "  " << curr->getInDegree();
//   // }
//   //
//   // std::cout << std::endl << std::endl;
//
//
//   std::cout << std::endl << std::endl<< std::endl<<std::endl;
//
//   std::vector<SetObj> bounds(NUM_DIRS, SetObj(NUM_BLOCKS));
//   SetObj cone(NUM_BLOCKS);
//
//   for (int i = -1;i<NUM_PERIODS;++i){
//
//     if (i > -1){
//       std::cout << "\t  *** period " << i << " boundaries ***";
//       getSingleInitBoundaries(sol,i,bounds);
//     }
//     else{
//       std::cout << "\n\n\t     *** no boundaries ***";
//     }
//
//     displayGraph(WIDTH, sol,bounds, cone);
//
//     std::cin.get();
//   }
//
//   sol.x = std::vector<int>{0,0,0,0,0,
//                            0,0,0,0,1,
//                            2,0,0,1,1,
//                            2,2,2,1,1,
//                            2,2,2,2,2};
//
//    NUM_BLOCKS = sol.x.size();
//    NUM_PERIODS = *max_element(sol.x.begin(), sol.x.end())+1;
//    WIDTH = (int)sqrt(NUM_BLOCKS);
//    graph = Graph();
//
//    genTestGraph(NUM_BLOCKS,WIDTH,graph);
//
//    model->graph = graph;
//
//    std::cout << std::endl << std::endl<< std::endl<<std::endl;
//
//    // reset boundaries
//    for (int i = 0; i < NUM_DIRS; ++i){
//      bounds[i].clear(NUM_BLOCKS);
//    }
//
//    const Node* curr = model->graph.getNode(24);
//    std::vector<Node*> connected = std::vector<Node*>();
//
//    std::cout << "24 connected to:" << std::flush;
//
//    // iterate over both directions and check connected blocks
//    for (int dir = 0; dir < NUM_DIRS; ++dir){
//      const int num_connected = curr->getConnected(dir, connected);
//
//      // iterate over all connected nodes
//      for (size_t node = 0; node < num_connected; ++node){
//        int adj_idx = connected[node]->getID();
//        std::cout << " " << adj_idx << std::flush;
//      }
//     }
//
//     std::cout << std::endl;
//
//    //  int swap_block = 17;
//    int period = 2;
//
//    getSingleInitBoundaries(sol,2,bounds);
//
//    for (int i = -1; i < 50; ++i){
//      const int dir = getDirection();
//      const int swap_block = bounds[dir].getRandomElement(rng);
//      if (swap_block == -1){
//        displayGraph(WIDTH, sol, bounds, cone);
//       //  std::cout << "period empty!" << std::endl;
//        continue;
//      }
//      if (i > -1){
//        cone.addElement(swap_block);
//        displayGraph(WIDTH, sol, bounds, cone);
//       //  std::cout << "choosing block" << std::endl;
//        cone.clear();
//        std::cin.get();
//        if (dir == BACKWARD){
//          sol.x[swap_block]++;
//        }
//        else {
//          sol.x[swap_block]--;
//        }
//        recomputeBlockBoundary(sol, dir, bounds, swap_block, period);
//      }
//      else{
//        std::cout << "\n\n\t     *** no boundaries ***";
//      }
//
//      displayGraph(WIDTH, sol, bounds, cone);
//
//     //  if (sol.x[swap_block] == period){
//     //    std::cout << swap_block << " added to forward\n" << std::endl;
//     //  }
//     //  else{
//     //    std::cout << swap_block << " added to backward\n" << std::endl;
//     //  }
//
//      std::cin.get();
//    }
// }
//
//
// void LocalSearch::runTests(Sol_Int &sol){
//   int MAX_DEPTH = 5;
//   std::vector<int> CONE_TIPS{100,45};
//
//   sol.x = std::vector<int>{1,0,0,0,0,0,0,0,0,0,0,
//                            1,1,0,0,0,0,0,0,0,0,0,
//                            1,1,1,3,0,0,0,0,0,0,2,
//                            1,1,3,3,3,0,0,0,0,2,2,
//                            1,3,3,3,3,3,0,2,2,2,2,
//                            5,3,3,3,3,3,4,2,2,2,2,
//                            5,5,3,3,3,4,4,4,2,2,2,
//                            5,5,5,5,4,4,4,4,4,4,2,
//                            5,5,5,5,5,4,4,4,4,4,4,
//                            5,5,5,5,5,5,4,4,4,4,4,
//                            5,5,5,5,5,5,5,4,4,4,4};
//
//   int NUM_BLOCKS = sol.x.size();
//   int NUM_PERIODS = *max_element(sol.x.begin(), sol.x.end())+1;
//   int WIDTH = (int)sqrt(NUM_BLOCKS);
//   Graph graph = Graph();
//
//   genTestGraph(NUM_BLOCKS,WIDTH,graph);
//
//   model->graph = graph;
//
//   // for (int b = 0; b < NUM_BLOCKS; ++b){
//   //   if (b % WIDTH == 0){
//   //     std::cout << std::endl << std::endl;
//   //   }
//   //
//   //   const Node* curr = graph.getNode(b);
//   //
//   //   std::cout << "  " << curr->getInDegree();
//   // }
//   //
//   // std::cout << std::endl << std::endl;
//
//
//   std::cout << std::endl << std::endl<< std::endl<<std::endl;
//
//   std::vector<SetObj> bounds(NUM_DIRS, SetObj(NUM_BLOCKS));
//   SetObj cone(NUM_BLOCKS);
//
//   int cone_depth = 0;
//
//   for (int i = -1;i<NUM_PERIODS;++i){
//
//     if (i > -1){
//       std::cout << "\t  *** period " << i << " boundaries ***";
//       getInitBoundaries(sol,i,bounds);
//     }
//     else{
//       std::cout << "\n\n\t     *** no boundaries ***";
//     }
//
//     displayGraph(WIDTH, sol,bounds, cone);
//
//     std::cin.get();
//   }
//
//   sol.x = std::vector<int>{0,0,0,0,0,0,0,0,0,0,0,0,
//                            0,0,0,0,0,0,0,0,0,0,0,0,
//                            0,0,0,0,0,0,0,0,0,0,0,0,
//                            0,0,0,0,0,0,0,0,0,0,0,0,
//                            0,0,0,0,0,1,1,0,0,0,0,0,
//                            0,0,0,1,1,1,1,1,0,0,0,0,
//                            0,0,1,1,1,1,1,1,1,1,1,0,
//                            1,1,1,1,1,1,1,1,1,1,1,1,
//                            1,1,1,1,1,1,1,1,1,1,1,1,
//                            1,1,1,1,1,1,1,1,1,1,1,1,
//                            1,1,1,1,1,1,1,1,1,1,1,1,
//                            1,1,1,1,1,1,1,1,1,1,1,1};
//
//
//    NUM_BLOCKS = sol.x.size();
//    NUM_PERIODS = *max_element(sol.x.begin(), sol.x.end())+1;
//    WIDTH = (int)sqrt(NUM_BLOCKS);
//    graph = Graph();
//
//    genTestGraph(NUM_BLOCKS,WIDTH,graph);
//
//    model->graph = graph;
//
//    std::cout << std::endl << std::endl<< std::endl<<std::endl;
//
//    // reset boundaries
//    for (int i = 0; i < NUM_DIRS; ++i){
//      bounds[i].clear(NUM_BLOCKS);
//    }
//
//    for (int i = 0;i<MAX_DEPTH+1;++i){
//
//      if (i > 0){
//        std::cout << "\t      *** depth " << i << " boundaries ***";
//        getInitBoundaries(sol,1,bounds);
//        increaseBoundaryDepth(sol,i,FORWARD,bounds[FORWARD]);
//        increaseBoundaryDepth(sol,i,BACKWARD,bounds[BACKWARD]);
//      }
//      else{
//        std::cout << "\n\n\t       *** no boundaries ***";
//      }
//
//      displayGraph(WIDTH, sol,bounds, cone);
//
//      std::cin.get();
//    }
//
//    std::cout << "\t      *** initial boundaries ***";
//    getInitBoundaries(sol,1,bounds);
//    displayGraph(WIDTH, sol,bounds, cone);
//
//    std::cin.get();
//
//    std::cout << "\t      *** depth " << 3 << " boundaries ***";
//    getInitBoundaries(sol,1,bounds);
//    increaseBoundaryDepth(sol,3,FORWARD,bounds[FORWARD]);
//    increaseBoundaryDepth(sol,3,BACKWARD,bounds[BACKWARD]);
//
//    displayGraph(WIDTH, sol,bounds, cone);
//
//    std::cin.get();
//
//    int dir = BACKWARD;
//
//    for (int i = 0; i < CONE_TIPS.size(); ++i){
//
//      std::cout << "      *** select random block in boundary ***";
//      cone.addElement(CONE_TIPS[i]);
//      displayGraph(WIDTH, sol,bounds, cone);
//      std::cin.get();
//
//      std::vector<double> cone_res_use(NUM_DIRS,0.0);
//      double cone_profit = 0.0;
//      int swap_period = -1;
//      cone.clear(NUM_BLOCKS);
//      // computeCone(sol, CONE_TIPS[i], cone, cone_profit, cone_res_use, swap_period, dir, cone_depth);
//      std::cout << "\t   *** compute cone of block (depth: " << cone_depth <<")***";
//      displayGraph(WIDTH, sol,bounds, cone);
//      std::cin.get();
//
//      std::cout << "\t      *** swap cone period ***";
//      // swapConePeriod(sol,swap_period,cone);
//      // clear cone data structure
//     //  cone.clear(NUM_BLOCKS);
//      displayGraph(WIDTH, sol,bounds, SetObj(NUM_BLOCKS));
//      std::cin.get();
//
//      recomputeConeBoundary(sol, 3, dir, bounds, cone);
//      displayGraph(WIDTH, sol,bounds, SetObj(NUM_BLOCKS));
//      std::cin.get();
//
//      cone.clear();
//
//     //  std::cout << "\t      *** show cone bounds ***";
//     //  getConeBoundaries(sol, cone, bounds,bounds);
//     //  // clear cone data structure
//     //  cone.clear();
//     //  displayGraph(WIDTH, sol,bounds, cone);
//     //  std::cin.get();
//      //
//     //  std::cout << "       *** recompute initial boundaries ***";
//     //  recomputeInitBoundary(sol,1,bounds);
//      //
//     //  displayGraph(WIDTH, sol,bounds, cone);
//     //  std::cin.get();
//      //
//     //  std::cout << "\t   *** new depth " << 3 << " boundaries ***";
//     //  getInitBoundaries(sol,1,bounds);
//     //  increaseBoundaryDepth(sol,3,FORWARD,bounds[FORWARD]);
//     //  increaseBoundaryDepth(sol,3,BACKWARD,bounds[BACKWARD]);
//      //
//     //  displayGraph(WIDTH, sol,bounds, cone);
//      //
//     //  std::cin.get();
//
//      dir = FORWARD;
//   }
//
//   sol.x = std::vector<int>{0,1,0,0,0,1,
//                            1,1,1,0,1,2,
//                            1,1,2,1,2,2,
//                            2,4,4,3,3,3,
//                            4,4,4,4,3,3,
//                            4,5,5,5,3,4};
//
//
//
//
//    // reset boundaries
//    for (int i = 0; i < NUM_DIRS; ++i){
//      bounds[i].clear(NUM_BLOCKS);
//    }
//    cone.clear(NUM_BLOCKS);
//
//    NUM_BLOCKS = sol.x.size();
//    NUM_PERIODS = *max_element(sol.x.begin(), sol.x.end())+1;
//    WIDTH = (int)sqrt(NUM_BLOCKS);
//    graph = Graph();
//
//    genTestGraph(NUM_BLOCKS,WIDTH,graph);
//
//    model->graph = graph;
//
//    displayGraph(WIDTH, sol,bounds, cone);
//
//    std::cin.get();
//
//    int cone_tip = 35;
//    std::vector<double> cone_res_use(NUM_DIRS,0.0);
//    double cone_profit = 0.0;
//    int swap_period = -1;
//    cone.clear(NUM_BLOCKS);
//    // computeCone(sol, cone_tip, cone, cone_profit, cone_res_use, swap_period, BACKWARD, cone_depth);
//
//    displayGraph(WIDTH, sol,bounds, cone);
//    std::cin.get();
//
//    // swapConePeriod(sol,swap_period,cone);
//
//    displayGraph(WIDTH, sol,bounds, cone);
//
//    std::cin.get();
//
//    cone_tip = 18;
//    swap_period = -1;
//    cone.clear(NUM_BLOCKS);
//    // computeCone(sol, cone_tip, cone, cone_profit, cone_res_use, swap_period, FORWARD,cone_depth);
//
//    swap_period = 4; // work out why that isnt working!
//
//    displayGraph(WIDTH, sol,bounds, cone);
//    std::cin.get();
//
//    // swapConePeriod(sol,swap_period,cone);
//
//    displayGraph(WIDTH, sol,bounds, cone);
//
//    std::cin.get();
//
//    cone_tip = 5;
//    swap_period = -1;
//    cone.clear(NUM_BLOCKS);
//    // computeCone(sol, cone_tip, cone, cone_profit, cone_res_use, swap_period, FORWARD,cone_depth);
//
//    swap_period = 2; // work out why that isnt working!
//
//    displayGraph(WIDTH, sol,bounds, cone);
//    std::cin.get();
//
//    // swapConePeriod(sol,swap_period,cone);
//
//    displayGraph(WIDTH, sol,bounds, cone);
//
//    std::cin.get();
//
//    cone_tip = 0;
//    swap_period = -1;
//    cone.clear(NUM_BLOCKS);
//    // computeCone(sol, cone_tip, cone, cone_profit, cone_res_use, swap_period, FORWARD,cone_depth);
//
//    swap_period = 1; // work out why that isnt working!
//
//    displayGraph(WIDTH, sol,bounds, cone);
//    std::cin.get();
//
//    // swapConePeriod(sol,swap_period,cone);
//
//    displayGraph(WIDTH, sol,bounds, cone);
//
//    std::cin.get();
//
//
//   // increaseBoundaryDepth(sol,TEST_PERIOD,BOUND_DEPTH,FORWARD,bound_bools[FORWARD],bound_sets[FORWARD]);
//   //
//   // displayGraph(WIDTH, sol,bound_bools);
//   //
//   // increaseBoundaryDepth(sol,TEST_PERIOD,BOUND_DEPTH,BACKWARD,bound_bools[BACKWARD],bound_sets[BACKWARD]);
//   //
//   // displayGraph(WIDTH, sol,bound_bools);
//
// }
//
//
// void LocalSearch::displayGraph(int width, Sol_Int &sol,
//   std::vector<SetObj> &bounds, const SetObj &cone){
//
//   //std::cout << "**** displaying graph ****\n" << std::endl;
//
//   const int nB = model->graph.getNumNodes();
//
//   for (int b = 0; b < nB; ++b){
//     if (b % width == 0){
//       std::cout << std::endl << "  ";
//       if (b > 0){
//         for (int i = b; i < b+width; ++i){
//           std::cout << "";
//           if (sol.x[i] != sol.x[i-width]){
//             std::cout << "\e[90m\x1b(0\x71\x1b(B\x1b(0\x71\x1b(B\x1b(0\x71\x1b(B\e[0m";
//           }
//           else{
//             std::cout << "   ";
//           }
//           std::cout << " ";
//         }
//       }
//       std::cout << std::endl;
//     }
//
//     std::cout << " ";
//
//     if (b % width != 0 && sol.x[b-1] != sol.x[b]){
//       std::cout << "\e[90m\x1b(0\x78\x1b(B\e[0m";
//     }
//     else{
//       std::cout << " ";
//     }
//
//     // if (b < 10)
//       std::cout << " ";
//
//
//     switch (sol.x[b]) {
//       case 0:{
//         if (bounds[0].is_element(b) || bounds[1].is_element(b)){
//           std::cout << "\e[30m\e[101m"; break;
//         }
//         else{std::cout << "\e[31m\e[1m"; break;}
//       }
//       case 1:{
//         if (bounds[0].is_element(b) || bounds[1].is_element(b)){
//           std::cout << "\e[30m\e[106m"; break;
//         }
//         else{std::cout << "\e[36m\e[1m"; break;}
//       }
//       case 2: {
//         if (bounds[0].is_element(b) || bounds[1].is_element(b)){
//           std::cout << "\e[30m\e[103m"; break;
//         }
//         else{std::cout << "\e[33m\e[1m"; break;}
//       }
//       case 3: {
//         if (bounds[0].is_element(b) || bounds[1].is_element(b)){
//           std::cout << "\e[30m\e[102m"; break;
//         }
//         else{std::cout << "\e[32m\e[1m"; break;}
//       }
//       case 4: {
//         if (bounds[0].is_element(b) || bounds[1].is_element(b)){
//           std::cout << "\e[30m\e[105m"; break;
//         }
//         else{std::cout << "\e[35m\e[1m"; break;}
//       }
//       case 5: {
//         if (bounds[0].is_element(b) || bounds[1].is_element(b)){
//           std::cout << "\e[30m\e[107m"; break;
//         }
//         else{std::cout << "\e[37m\e[1m"; break;}
//       }
//     }
//
//     if (cone.is_element(b)){
//       if (sol.x[b] != 3){
//         std::cout << "\e[0m\e[30m\e[102m";
//       }
//       else{
//         std::cout << "\e[0m\e[30m\e[101m";
//       }
//     }
//
//     // std::cout << b << "\e[0m";
//
//     std::cout << sol.x[b] << "\e[0m";
//
//   }
//   std::cout << std::endl << std::endl<< std::endl << std::endl;
// }
//
// void LocalSearch::genTestGraph(int nB, int width, Graph &graph){
//
//
//   if (std::floor((float)nB/width) != (float)nB/width){
//     std::cout << "nB/height not integer!" << std::endl;
//     return;
//   }
//
//   // adding nodes
//   for (int b = 0; b < nB; ++b){
//     graph.addNode(b);
//   }
//
//   // adding arcs
//   int arcID = 0;
//   for (int y = 0; y < (nB/width)-1; y++){
//     for (int x = 0; x < width; ++x){
//       if (x > 0){
//         graph.addArc(arcID++, (y*width)+x, ((y+1)*width)+(x-1));
//       }
//       graph.addArc(arcID++, (y*width)+x, ((y+1)*width)+x);
//       if (x < width-1){
//         graph.addArc(arcID++, (y*width)+x, ((y+1)*width)+(x+1));
//       }
//     }
//   }
// }

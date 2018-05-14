#include "RandomSearch.h"

double RandomSearch::betterRandomSearch(const BranchNode_info &probInfo, Sol_Int &sol,SinglePModel *model){
  const int nB = model->getNBlock();
  const int t_max = model->getNPeriod();
  const int d_max = model->getnDestination();
  const int r_max = model->getnResources();
  const double rate = model->getDiscountRate();

  std::vector<Block> * blocks=model->getBlock();

  std::cout << "\n\nnumber of blocks in problem: " << nB << std::endl << std::endl;

  //std::cout <<"\n***** RUNNING RANDOM SEARCH *****\n" << std::endl;

  // std::random_device rd;     // only used once to initialise (seed) engine
  // std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)

  int block_cnt = 0;

  long double mine_total = 0.0;

  //std::vector<bool> open_blocks(nB, false);
  std::vector<int> open_idx;
  std::vector<bool> mined(nB, false);

  for (size_t b=0;b<nB;++b){
    if ((*blocks)[b].getNumPred() == 0 && probInfo.time[b][0] < t_max){
      //open_blocks[b] = true;
      open_idx.push_back(b);
      //std::cout << b << std::endl;
    }
  }

  // check initial list for invalid blocks
  for (int i = 0; i< open_idx.size(); ++i){
    if (open_idx[i] >= nB) {
      std::cout << "******** block " << open_idx[i] << " detected initially!! ***********" << std::endl;
    }
  }

  bool stop_mining = false;

  int period = 0;

  std::vector<double> res_limit (r_max, 0.0);
  std::vector<long double> res_use (r_max, 0.0);
  long double p_res_use = 0.0;

  // compute resource used for mining blocks
  int r_mining = 0;
  for (size_t r = 0;r < r_max; ++r){
    res_limit[r] = model->getLimit(r, 0);
  }

  bool run_tests = false;

  //stop_mining = true;

  // multimap to store all blocks mined in period, sorted by profit difference
  std::multimap<double, int> p_mined;
  std::multimap<double, int>::reverse_iterator i;

  while (!stop_mining && open_idx.size()>0){
    std::uniform_int_distribution<int> uni(0,open_idx.size()-1); // guaranteed unbiased
    int curr_idx = uni(rng); // index for open_idx

    // get block index associated with current open_idx index
    int b = open_idx[curr_idx];

    // if (b >= nB){ std::cout << "**** DANGER!!! block " << b  << " detected!! ****" << std::endl;
    //
    // for (int ii = 0;ii < open_idx.size(); ++ii) std::cout << open_idx[ii] << std::endl;
    // }
    // get block
    const Block & block=model->getBlock(b);
    //std::cout <<"testing block " << b << ", profit: " << std::endl;// << block.getProfit(0) << "..." <<std::flush;

    bool res_exceeded = false;

    // if any of the resources are exceeded (allowing for scaled resource use) do period mining
    for (int r = 0; r < r_max; ++r){
      if ((res_use[r] + block.getRCoef(0,r)) > res_limit[r]*sh.SCALED_RESOURCE){
        res_exceeded = true;
        // reset res_use vector
        std::fill(res_use.begin(), res_use.end(),0.0);
        break;
      }
    }

    // if adding block to list will exceed limit for current period, do period mining
    if (res_exceeded){
      // resource use vector for period mining
      std::vector<long double> p_res_use (r_max, 0.0);

      // iterate over multimap to process highest value blocks
  //      std::cout << "*** period " << period << " ***\nblocks processed: "<<std::flush;
      int p_count = 0;
      int t_count = 0;
      for(i = p_mined.rbegin();i != p_mined.rend();++i){
        const Block & c_block = model->getBlock((*i).second);

        bool p_res_exceeded = false;
        for (int r = 0; r < r_max; ++r){
          if ((p_res_use[r] += c_block.getRCoef(0,r)) > res_limit[r]){
            p_res_exceeded = true;
            break;
          }
        }
        if (p_res_exceeded){
          break;
        }

        double block_profit = c_block.getProfit(0) / pow(1+rate,period);
        mine_total += block_profit;
        sol.x[(*i).second] = period;
        p_count++;
        t_count++;
      }

      //start new period
      period++;

      if (period > t_max){
        stop_mining=true;
        break;
      }

      p_mined.clear();
    }
    // *** end period mining ***

    if (stop_mining){
      break;
    }

    // add block to current periods mining list
    p_mined.insert(p_mined_pair((block.getProfit(0)), b));

    for (int r = 0; r < r_max; ++r){
      res_use[r] += block.getRCoef(0,r);
    }

    if (curr_idx < open_idx.size()-1){
      open_idx[curr_idx] = open_idx.back();
    }
    open_idx.pop_back();

    // mark block as mined
    mined[b] = true;
    block_cnt++;

    const Node* curr = model->graph.getNode(b);
    for (size_t next_arc = 0; next_arc < curr->getOutDegree(); ++next_arc){
      int next_idx = curr->getOutArc(next_arc)->getTgtID();
      if (next_idx >= nB) continue;

      const Node* next = model->graph.getNode(next_idx);
      bool block_open = true;
      for (size_t prev_arc = 0; prev_arc < next->getInDegree(); ++prev_arc){
        int prev_idx = next->getInArc(prev_arc)->getSrcID();
        if (!mined[prev_idx]){
          block_open = false;
          break;
        }
      }
      if (block_open && !mined[next_idx] && probInfo.time[next_idx][0] < t_max){
        open_idx.push_back(next_idx);
        //if (next_idx == 1060){ std::cout << "**** DANGER!!! ADDING 1060 ****" << std::endl;}
      }
    }
  }

  std::cout << "\nbetterRandom search completed!\n\n" << block_cnt << " blocks mined\n"
            << mine_total << " total mine value\n\n";

  //std::cout << "blocks not mined: " << std::flush;
  //for (int i = 0; i < mined.size(); ++i){
  //  if (!mined[i]) std::cout << i << " ";
  //}
  //std::cout << std::endl;


  sol.nT = t_max;
  sol.obj = mine_total;
  // try{
  //   std::cout << "\nChecking solution from solver...\n" << std::endl;
  //
  //   bool test_error = verify((*model), sol);
  //
  //   std::cout << "Solution found by solver was ";
  //   if (!test_error)
  //     std::cout << "feasible!\n" << std::endl;
  //   else{
  //     std::cout << "infeasible!\n" << std::endl;
  //     throw qol::Exception("Infeasible solution!");
  //   }
  // } // end try statement
  // catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
  // catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

 return mine_total;
}

double RandomSearch::randomSearch(const BranchNode_info &probInfo, Sol_Int &sol, SinglePModel *model){
  const int nB = model->getNBlock();
  const int t_max = model->getNPeriod();
  const int d_max = model->getnDestination();
  const int r_max = model->getnResources();
  const double rate = model->getDiscountRate();

  std::vector<Block> * blocks=model->getBlock();

  //std::random_device rd;     // only used once to initialise (seed) engine
  //std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)

  int block_cnt = 0;

  long double mine_total = 0.0;

  std::vector<int> open_idx;
  std::vector<bool> mined(nB, false);

  for (size_t b=0;b<nB;++b){
    if ((*blocks)[b].getNumPred() == 0 && probInfo.time[b][0] < t_max){
      open_idx.push_back(b);
    }
  }

  bool stop_mining = false;

  int period = 0;

  std::vector<double> res_limit (r_max, 0.0);
  std::vector<long double> res_use (r_max, 0.0);

  for (size_t r = 0;r < r_max; ++r){
    res_limit[r] = model->getLimit(r, 0);
  }

  while (!stop_mining && open_idx.size()>0){
    std::uniform_int_distribution<int> uni(0,open_idx.size()-1); // guaranteed unbiased
    int curr_idx = uni(rng); // index for open_idx

    // get block index associated with current open_idx index
    int b = open_idx[curr_idx];

    // if block is being mined too early, leave it and go to next
    if (probInfo.time[b][1] < period){
      std::cout << "too early!" << std::endl;
      continue;
    }

    // if block is being mined too late, remove from list and go to next
    if (probInfo.time[b][0] > period){
      std::cout << "too late!" << std::endl;
      // remove from open_idx list by swapping with end element
      if (curr_idx < open_idx.size()-1){
        open_idx[curr_idx] = open_idx.back();
      }
      open_idx.pop_back();
      continue;
    }

    // get block
    const Block & block=model->getBlock(b);

    double block_profit = block.getProfit(0) / pow(1+rate,period);

    bool change_period = false;
    // increase resource usage and check if resource limit reached
    for (int r = 0; r < r_max; ++r){
      if ((res_use[r] + block.getRCoef(0,r)) > res_limit[r] * sh.SCALED_RESOURCE){
        // store res_use vector and reset, then increment period
        sol.res_use[period] = res_use;
        for (int r2 = 0; r2 < r_max; ++r2){
          res_use[r2] = block.getRCoef(0,r2);
        }
        period++;
        if (period >= t_max){
          stop_mining = true;
        }
        break;
      }
      else{
        res_use[r] += block.getRCoef(0,r);
      }
    }
    // if resource limit reached, change period
    if (stop_mining){
      break;
    }
    // increment block counter
    block_cnt++;

    // increase total profit
    mine_total += block.getProfit(0)/pow(1+rate,period);

    // remove from open_idx list by swapping with end element
    if (curr_idx < open_idx.size()-1){
      open_idx[curr_idx] = open_idx.back();
    }
    open_idx.pop_back();

    // mark block as mined
    mined[b] = true;

    sol.x[b] = period;

    const Node* curr = model->graph.getNode(b);
    for (size_t next_arc = 0; next_arc < curr->getOutDegree(); ++next_arc){
      int next_idx = curr->getOutArc(next_arc)->getTgtID();
      if (next_idx >= nB) continue;

      const Node* next = model->graph.getNode(next_idx);
      bool block_open = true;
      for (size_t prev_arc = 0; prev_arc < next->getInDegree(); ++prev_arc){
        int prev_idx = next->getInArc(prev_arc)->getSrcID();
        if (!mined[prev_idx]){
          block_open = false;
          break;
        }
      }
      if (block_open && !mined[next_idx] && probInfo.time[next_idx][0] < t_max){
        open_idx.push_back(next_idx);
      }
    }
  }

  std::cout << "\nRandom search completed!\n\n" << block_cnt << " blocks mined\n"
            << mine_total << " total mine value\n\n";

  sol.nT = t_max;
  sol.obj = mine_total;
  // try{
  //   std::cout << "\nChecking solution from solver...\n" << std::endl;
  //
  //   bool test_error = verify((*model), sol);
  //
  //   std::cout << "Solution found by solver was ";
  //   if (!test_error)
  //     std::cout << "feasible!\n" << std::endl;
  //   else{
  //     std::cout << "infeasible!\n" << std::endl;
  //     throw qol::Exception("Infeasible solution!");
  //   }
  // } // end try statement
  // catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
  // catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}

 return mine_total;
}

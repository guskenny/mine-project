#include "ConeMiner.h"
#include "UpitSolver.h"
#include <algorithm>
#include <atomic>

void ConeMiner::solve(Sol_Int &sol){

  const int t_max = prob.getNPeriod();
  const int r_max = prob.getnResources();
  const double rate = prob.getDiscountRate();
  const int nB = prob.getNBlock();

  std::vector<bool> mined(nB,0);
  double mine_total = 0.0;
  int mined_cnt = 0;

  for (int t = 0; t < t_max; ++t){
    // std::cout << "period: " << t << std::endl;
    int period_cnt = 0;
    // std::cout << "fixing earliest... " << std::flush;
    fixEarliest(mined, t);
    // std::cout << "done!" << std::endl << "getting cones... " << std::flush;
    std::vector<bool> included(nB,0);
    std::vector<double> res_limit (r_max,0.0);
    std::vector<double> res_use (r_max,0.0);
    for (int i = 0; i < sh.NUM_SUB_PERIODS; ++i){
      std::cout << "*** Subperiod " << i << " of " << sh.NUM_SUB_PERIODS << " ***" << std::endl;
      getMostValuableCones(mined, included, res_limit, res_use, mine_total, t);
    }
    // std::cout << "done!" << std::endl << "finishing randomly... " << std::flush;
    finishRandom(mined, included, res_limit, res_use, mine_total, t);
    // std::cout << "done!" << std::endl << "updating solution... " << std::flush;
    for (int b = 0; b < nB; ++b){
      if (included[b]){
        sol.x[b] = t;
        mined_cnt++;
        period_cnt++;
      }
      // std::cout << "done!" << std::endl;
    }

    std::cout << "period " << t << ": " << period_cnt << " blocks mined (" << mined_cnt << " total)" << std::endl;
  }

  sol.nT = t_max;
  sol.obj = mine_total;

  std::cout << "Cone mining complete! " << mined_cnt << " blocks mined for total value of " << mine_total << std::endl;

}

// the earliest a block can be mined is if we have enough resources to clear
// all of the cone of predecessor blocks by that period
int ConeMiner::fixEarliest(const std::vector<bool> &mined, const int period)
{
    //std::cout << "\tPreprocess::fixEarliest disabled!!!!\n"; return;
    const double eps=1e-6;
    const int t_max = prob.getNPeriod();
    const int d_max = prob.getnDestination();
    const int r_max = prob.getnResources();
    const double rate = prob.getDiscountRate();
    const int nB = prob.getNBlock();

    int cont=0;

    // for (int i=0;i<nB;++i)
    //   if (mined[i])
    //     cont++;
    //
    // std::cout << cont << " total blocks mined\n"<<std::endl;

    cumRes.clear();
    cumRes.resize(t_max+1,std::vector<long double> (r_max,0.0));

    // set up vector of cumulative resource limits for each time period and each resource
    for(int t=0;t<t_max;++t){
	for(int r=0;r<r_max;++r){
	    cumRes[t][r] += prob.getLimit(r,t);
	    cumRes[t+1][r]= cumRes[t][r];

	}
    }

    std::atomic<int> cnt; cnt=0;
    bool infeas=false;

    // establish vectors for cone values, cone membership and resource use for each block
    cone_value.clear();
    cone_value.resize(nB, 0.0);
    inCone.clear();
    inCone.resize(nB, std::vector<bool> (nB, false)); // true if in cone
    res.clear();
    res.resize(nB, std::vector<long double> (prob.getnResources(), 0.0));

    double max_cone_value = -9e10;

//#   pragma omp parallel for
    for(int b=0;b<nB;++b){
	std::vector<int> stack; // stack for search
	stack.push_back(b);
	inCone[b][b]=true; // set current block to in cone
	while(! stack.empty()){
	    int blockID = stack.back(); // get next block id
	    const Block & block=prob.getBlock(blockID); // get block
	    stack.pop_back(); // remove block from stack

      //double max_profit = block.getProfit(0);
      double max_profit = block.getProfit(0)/pow(1+rate,period);
      // for (int d=1;d<d_max;++d){
      //   //max_profit = std::max(max_profit, block.getProfit(d));
      //   max_profit = std::max(max_profit, block.getProfit(d)/pow(1+rate,period));
      // }

      cone_value[b] += max_profit;
	    for(int r=0;r<r_max;++r){ // iterate over resources
		long double minR=block.getRCoef(0,r); // get resource from destination 0
		// for(int d=1;d<d_max;++d)
		//     minR = std::min(minR,(long double)block.getRCoef(d,r)); // find minimum resource
		res[b][r] += minR; // add minimum resource use to resource vector
	    }
      // iterate over predecessors
	    for(auto p=block.getPreds().begin();p!=block.getPreds().end();++p)
        if( ! inCone[b][*p] && !mined[*p] && (period > node.fixed[*p][0])){ // if not in cone and not previously mined
		      inCone[b][*p] = true; // add to cone
		      stack.push_back(*p); // push back onto stack
		}
	}

  if (cone_value[b] > max_cone_value)
    max_cone_value = cone_value[b];
	int earliest=node.time[b][0]; // get previously defined earliest mining time
  //if (earliest < 0) earliest=0; // REMOVE EVENTUALLY!!!!
  //std::cout << earliest << " " << std::endl;
	for(int r=0;r<r_max;++r){
      // while resource used is greater than cumulative available, increased earliest
	    while(earliest < t_max && res[b][r] > (1+eps)*cumRes[earliest][r]){
		      ++earliest;
      }
  }
	cnt += earliest - node.time[b][0]; // count how many time periods removed
	node.time[b][0] = earliest; // update value in BranchNode_info
	//check if earliest mineable time is greater than latest mineable time
  if(earliest >= node.time[b][1] && node.time[b][1] < t_max)
	  infeas=true;
    }



    std::cout << "fixEarliest removed " << cnt << " periods, avg "
	      << double(cnt)/nB << " per block\n";
    if(infeas) cnt = Infeasible;
    return cnt;
}

// the earliest a block can be mined is if we have enough resources to clear
// all of the cone of predecessor blocks by that period
void ConeMiner::getMostValuableCones(std::vector<bool> &mined, std::vector<bool> &included, std::vector<double> &res_limit, std::vector<double> &res_use, double &mine_total, const int t)
{
    //std::cout << "\tPreprocess::fixEarliest disabled!!!!\n"; return;
    //const double eps=1e-6;
    //const int t_max = prob.getNPeriod();
    const int r_max = prob.getnResources();
    const int nB = prob.getNBlock();
    const double rate = prob.getDiscountRate();

    std::vector<double> res_to_use (r_max,0.0);

    int num_top_list = nB * sh.TOP_PERCENT;

    std::priority_queue<std::pair<double, int>> q;

    double list_min = 9e10;

    for (size_t r=0;r<r_max;++r){
      res_limit[r] = prob.getLimit(r,t) * double(1.0/sh.NUM_SUB_PERIODS);
    }

    // populate the initial num_top_list entries then only if bigger than smallest on list
    for (int b = 0; b<nB; ++b){
      bool cutoff_reached = false;
      for (size_t r = 0;r<r_max;++r){
        if (res[b][r] < res_limit[r] * sh.RES_CUTOFF)
          cutoff_reached = true;
      }
      if ((cone_value[b] > list_min || b < num_top_list)
              && t >= node.time[b][0] && t < node.time[b][1]
              && cone_value[b] > 0 && cutoff_reached){
        q.push(std::pair<double, int>(cone_value[b],b));
        if (cone_value[b] < list_min)
          list_min = cone_value[b];
      }
    }

    std::vector<int> top_list;

    top_list.reserve(num_top_list);

    top_list.push_back(q.top().second);
    q.pop();

    std::cout << "top list size: " << top_list.size() << std::endl;

    while (q.size() > 0){
      bool skip = false;
      for (int i = 0;i<top_list.size();++i){
        if (inCone[top_list[i]][q.top().second]){
          skip = true;
          break;
        }
      }
      if (!skip)
        top_list.push_back(q.top().second);
      q.pop();
    }


    std::cout << top_list.size() << " positive independent cones available in period "
              << t << std::endl;

    bool over_limit = false;

    int block_count = 0;
    int cone_count = 0;

    std::vector<bool> selected(top_list.size(),false);

    int cone = 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    // iterate over top valued cones until resource limit reached
    //std::cout << "selecting cones.. " << std::flush;
    while (!over_limit && cone_count < top_list.size()){
      std::vector<bool> to_include(nB, false);
      //std::cout << "selecting cone.. " << std::flush;
      if (dis(gen) <= sh.P_SELECTION && !selected[cone]){
        //std::cout << "cone " << cone << " selected!" << std::endl;
        cone_count++;
        //std::cout << cone_count << " cones selected!" <<std::endl;
        selected[cone] = true;
	      int root_b = top_list[cone];
        //std::cout << "root_b: " << root_b << std::endl;
	      for (size_t b=0;b<nB;++b){
    //        std::cout << "blerp.." << inCone[root_b][b] << std::flush;
	        if (inCone[root_b][b] && !to_include[b] && !mined[b]){
            //std::cout << "yup " << std::flush;
	          block_count++;
	          to_include[b]=true;
		        const Block & block=prob.getBlock(b); // get block
            if (r_max > 1){
	    	      for(int r=0;r<r_max;++r){ // iterate over resources
			          long double minR=block.getRCoef(0,r); // get resource from destination 0
			          // for(int d=1;d<d_max;++d){
	              //   // find minimum resource
			          //   minR = std::min(minR,(long double)block.getRCoef(d,r));
	              // }
			          res_to_use[r] += minR; // add minimum resource use to resource vector
		          }
            }
            else{
              res_to_use[0] += block.getRCoef(0,0);
            }
	        }
      //      std::cout << " blorp" << std::endl;
	      }
	      for (size_t r=0;r<r_max;++r){
          //std::cout << "resource " << r << ": " << res_to_use[r]
          //          << " of " << res_limit[r] << std::endl;
	        if (res_to_use[r] > (res_limit[r] * sh.LIMIT_MULTIPLIER))
	          over_limit = true;
	      }
        if (!over_limit){
          // include in the list of blocks
          for (int b = 0; b < nB; ++b){
            if (to_include[b]){
              mined[b] = true;
              included[b] = true;
              const Block & block=prob.getBlock(b);
              mine_total += block.getProfit(0) / pow(1+rate,t);
            }
          }
          for (size_t r=0;r<r_max;++r){
            res_use[r] = res_to_use[r];
          }
        }
			}
      //else{
        //std::cout << "cone " << cone << " skipped!" << std::endl;
      //}
      cone++;
			if (cone >= top_list.size()){
        cone = 0;
      }
    }
    //std::cout << "done!" << std::endl;

//    // iterate over top valued cones until resource limit reached
//    for (size_t cone=0;cone<top_list.size() && !over_limit;++cone){
//      cone_count++;
//      int root_b = top_list[cone];
//      for (size_t b=0;b<nB;++b){
//        if (inCone[root_b][b] && !included[b]){
//          block_count++;
//          included[b]=true;
//	        const Block & block=prob.getBlock(b); // get block
//    	    for(int r=0;r<r_max;++r){ // iterate over resources
//		        long double minR=block.getRCoef(0,r); // get resource from destination 0
//		        for(int d=1;d<d_max;++d){
//              // find minimum resource
//		          minR = std::min(minR,(long double)block.getRCoef(d,r));
//            }
//		        res_use[r] += minR; // add minimum resource use to resource vector
//	        }
//        }
//      }
//      for (size_t r=0;r<r_max;++r){
//        if (res_use[r] > (res_limit[r] * sh.LIMIT_MULTIPLIER))
//          over_limit = true;
//      }
//    }

    std::cout << block_count << " blocks added from " << cone_count << " cones" << std::endl;
}

void ConeMiner::finishRandom(std::vector<bool> &mined, std::vector<bool> &included, std::vector<double> &res_limit, std::vector<double> &res_use, double &mine_total, const int period){

  const int nB = prob.getNBlock();
  const int t_max = prob.getNPeriod();
  const int d_max = prob.getnDestination();
  const int r_max = prob.getnResources();
  const double rate = prob.getDiscountRate();

  //std::random_device rd;     // only used once to initialise (seed) engine
  //std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)

  int block_cnt = 0;

  std::vector<int> open_idx;
  // std::vector<bool> mined(nB, false);

  // find open blocks
  for (size_t b=0;b<nB;++b){
    if (mined[b]){
      continue;
    }
    const Block & block=prob.getBlock(b);
    if (block.getNumPred() == 0 && node.time[b][0] < t_max){
      open_idx.push_back(b);
      continue;
    }

    bool not_mined = false;
    for (auto p = block.getPreds().begin(); p != block.getPreds().end(); ++p){
      if (!mined[*p]){
        not_mined = true;
        break;
      }
    }
    if (!not_mined){
      open_idx.push_back(b);
    }
  }

  bool stop_mining = false;

  while (!stop_mining && open_idx.size()>0){
    std::uniform_int_distribution<int> uni(0,open_idx.size()-1); // guaranteed unbiased
    int curr_idx = uni(rng); // index for open_idx

    // get block index associated with current open_idx index
    int b = open_idx[curr_idx];

    // // if block is being mined too early, leave it and go to next
    // if (probInfo.time[b][1] < period){
    //   std::cout << "too early!" << std::endl;
    //   continue;
    // }

    // if block is being mined too late, remove from list and go to next
    if (node.time[b][0] > period || node.time[b][1] < period){
      // std::cout << "too late or early!!" << std::endl;
      // remove from open_idx list by swapping with end element
      if (curr_idx < open_idx.size()-1){
        open_idx[curr_idx] = open_idx.back();
      }
      open_idx.pop_back();
      continue;
    }

    // get block
    const Block & block=prob.getBlock(b);

    double block_profit = block.getProfit(0) / pow(1+rate,period);

    bool change_period = false;
    // increase resource usage and check if resource limit reached
    for (int r = 0; r < r_max; ++r){
      if ((res_use[r] + block.getRCoef(0,r)) > res_limit[r] * sh.SCALED_RESOURCE){
        stop_mining = true;
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
    included[b] = true;

    const Node* curr = prob.graph.getNode(b);
    for (size_t next_arc = 0; next_arc < curr->getOutDegree(); ++next_arc){
      int next_idx = curr->getOutArc(next_arc)->getTgtID();
      if (next_idx >= nB) continue;

      const Node* next = prob.graph.getNode(next_idx);
      bool block_open = true;
      for (size_t prev_arc = 0; prev_arc < next->getInDegree(); ++prev_arc){
        int prev_idx = next->getInArc(prev_arc)->getSrcID();
        if (!mined[prev_idx]){
          block_open = false;
          break;
        }
      }
      if (block_open && !mined[next_idx] && node.time[next_idx][0] < t_max){
        open_idx.push_back(next_idx);
      }
    }
  }

  //
  std::cout << block_cnt << " random blocks mined\n" << std::endl;
  // sol.nT = t_max;
  // sol.obj = mine_total;
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
}

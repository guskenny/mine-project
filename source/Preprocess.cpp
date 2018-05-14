#include "Preprocess.h"
#include "UpitSolver.h"
#include <algorithm>
#include <atomic>

void Preprocess::fixUPIT()
{
    UpitSolver upit(prob);
    upit.solve();
    const int t_max = prob.getNPeriod();
    std::vector<int> include(prob.getNBlock(),0);
    upit.getClosure(include);
    int cnt=0;
    for(int b=0;b<prob.getNBlock();++b)
		if(!include[b]){
			node.time[b][0]=node.time[b][1]=t_max;
			if(++cnt < 10) std::cout << "\t\tBlock " << b << " is never mined\n";
			node.mine_block[b] = false;
		}
    //residualProfit.resize(prob.getNBlock());

    //upit.getResidualProfit(residualProfit);
//    std::cout << "\tPreprocess::fixUPIT fixed " << cnt << " blocks, "
//			  << "profit upper bound = " << upit.getObjective(include) << std::endl;
}

int Preprocess::getProcessable(std::vector<bool> &processable){
    const int nB = prob.getNBlock();
    const int d_max = prob.getnDestination();
    const double rate = prob.getDiscountRate();

    std::vector<double> profits (nB,0.0);

    int b_count = 0;

    processable.clear();
    processable.resize(nB, true);

    // code only works for two destinations
    if (d_max != 2){
      return;
    }

    double max_profit = 0.0;

    for (size_t b=0;b<nB;++b){
	    const Block & block=prob.getBlock(b); // get block
      double process_profit = block.getProfit(1);
      profits[b] = process_profit;
      if (process_profit > max_profit)
        max_profit = process_profit;
    }

    for (size_t b=0;b<nB;++b){
      if (profits[b] < max_profit * sh.PROCESS_THRESHOLD){
        b_count++;
        processable[b] = false;
      }
    }


//    for (size_t b=0;b<nB;++b){
//      double no_process_profit = block.getProfit(0)/pow(1+rate,period);
//      double process_profit = block.getProfit(1)/pow(1+rate,period);
//      if ((process_profit-no_process_profit) < sh.PROCESS_THRESHOLD
//          || process_profit < 0){
//       processable[b] = false;
//       b_count++;
//      }
//    }

    return b_count;
}


// the earliest a block can be mined is if we have enough resources to clear
// all of the cone of predecessor blocks by that period
int Preprocess::fixEarliest(const std::vector<bool> &mined, const int period)
{
    //std::cout << "\tPreprocess::fixEarliest disabled!!!!\n"; return;
    const double eps=1e-6;
    const int t_max = prob.getNPeriod();
    const int d_max = prob.getnDestination();
    const int r_max = prob.getnResources();
    const double rate = prob.getDiscountRate();
    const int nB = prob.getNBlock();

    int cont=0;

    for (int i=0;i<nB;++i)
      if (mined[i])
        cont++;

    std::cout << cont << " total blocks mined\n"<<std::endl;

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
void Preprocess::getMostValuableCones(int t, std::vector<bool> &included, const std::vector<bool> &mined)
{
    //std::cout << "\tPreprocess::fixEarliest disabled!!!!\n"; return;
    //const double eps=1e-6;
    //const int t_max = prob.getNPeriod();
    const int d_max = prob.getnDestination();
    const int r_max = prob.getnResources();
    const int nB = prob.getNBlock();

    included.clear();
    included.resize(nB,false);


    int num_top_list = nB * sh.TOP_PERCENT;

    std::priority_queue<std::pair<double, int>> q;

    double list_min = 9e10;

    std::vector<double> res_limit (r_max,0.0);

    std::vector<double> res_use (r_max,0.0);

    for (size_t r=0;r<r_max;++r){
      res_limit[r] = prob.getLimit(r,t);
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
	        if (inCone[root_b][b] && !included[b]){
            //std::cout << "yup " << std::flush;
	          block_count++;
	          included[b]=true;
		        const Block & block=prob.getBlock(b); // get block
            if (r_max > 1){
	    	      for(int r=0;r<r_max;++r){ // iterate over resources
			          long double minR=block.getRCoef(0,r); // get resource from destination 0
			          // for(int d=1;d<d_max;++d){
	              //   // find minimum resource
			          //   minR = std::min(minR,(long double)block.getRCoef(d,r));
	              // }
			          res_use[r] += minR; // add minimum resource use to resource vector
		          }
            }
            else{
              res_use[0] += block.getRCoef(1,0);
            }
	        }
      //      std::cout << " blorp" << std::endl;
	      }
	      for (size_t r=0;r<r_max;++r){
          //std::cout << "resource " << r << ": " << res_use[r]
          //          << " of " << res_limit[r] << std::endl;
	        if (res_use[r] > (res_limit[r] * sh.LIMIT_MULTIPLIER))
	          over_limit = true;
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

    std::cout << block_count << " blocks added from " << cone_count << " cones"<<std::endl;
}


void Preprocess::fixLatest()
{
}

#ifndef __CONEMINER_H__
#define __CONEMINER_H__
//#include "MineProblem.h"
#include "SinglePSolver.h"
#include "SinglePModel.h"
#include "SettingsHandler.h"
#include <limits>
#include <queue>
#include <random>


// Preprocessor class tries to reduce set of possible times
// when a block might be mined
// this class modifies the BranchNode_info passed to it
// in by setting time[b][0/1] to earliest/latest time for mining b
class ConeMiner {
public:
    static constexpr int Infeasible = std::numeric_limits<int>::max();
    const SinglePModel &prob;
    BranchNode_info &node;
    const SettingsHandler &sh;

    std::mt19937 rng;

    std::vector<std::vector<long double> > cumRes;
    std::vector<double> cone_value;
    std::vector<std::vector<bool> > inCone;
    std::vector<std::vector<long double> > res;
    std::vector<int> skip_list;

    // by default Preprocess runs all fixing methods, this can be turned off by setting doIt=false
    ConeMiner(const SinglePModel &_prob, BranchNode_info &_node, const SettingsHandler &_sh) : prob(_prob), node(_node), sh(_sh) {
	    if( node.time.empty() ){ // not initialised
	       node.init(_prob.getNBlock(),_prob.getNPeriod(),prob.getnDestination(),
		       prob.getnResources());
      }
      std::random_device r;
      std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
      rng = std::mt19937(seed);
    }

    void solve(Sol_Int &sol);
    // earliest resource-feasible solution: returns no. fixed or INFEASIBLE
    int fixEarliest(const std::vector<bool> &mined, const int period);
    void getMostValuableCones(std::vector<bool> &mined, std::vector<bool> &included, std::vector<double> &res_limit, std::vector<double> &res_use, double &mine_total, const int period);
    void finishRandom(std::vector<bool> &mined, std::vector<bool> &included, std::vector<double> &res_limit, std::vector<double> &res_use, double &mine_total, const int period);
};


#endif

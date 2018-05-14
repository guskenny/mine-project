#ifndef __PREPROCESS_H__
#define __PREPROCESS_H__
//#include "MineProblem.h"
#include "SinglePSolver.h"
#include "SettingsHandler.h"
#include <limits>
#include <queue>


// Preprocessor class tries to reduce set of possible times
// when a block might be mined
// this class modifies the BranchNode_info passed to it
// in by setting time[b][0/1] to earliest/latest time for mining b
class Preprocess {
public:
    std::vector<double> residualProfit;
    static constexpr int Infeasible = std::numeric_limits<int>::max();
    const Daten &prob;
    BranchNode_info &node;
    const SettingsHandler &sh;

    std::vector<std::vector<long double> > cumRes;
    std::vector<double> cone_value;
    std::vector<std::vector<bool> > inCone;
    std::vector<std::vector<long double> > res;

    // by default Preprocess runs all fixing methods, this can be turned off by setting doIt=false
    Preprocess(const Daten &_prob, BranchNode_info &_node,const std::vector<bool> &mined, const SettingsHandler &_sh, bool doIt=true) : prob(_prob), node(_node), sh(_sh) {
	if( node.time.empty() ) // not initialised
	    node.init(_prob.getNBlock(),_prob.getNPeriod(),prob.getnDestination(),
		      prob.getnResources());
	if(doIt){
	  fixUPIT();
	  fixEarliest(mined,0);
	  fixLatest();
	}
    }
    void fixUPIT(); // eliminate blocks that are not in the UPIT solution
    // earliest resource-feasible solution: returns no. fixed or INFEASIBLE
    int fixEarliest(const std::vector<bool> &mined, const int period);
    int getProcessable(std::vector<bool> &processable);
    void fixLatest();	// latest completion time to do whole UPIT
    void getMostValuableCones(int t, std::vector<bool> &included, const std::vector<bool> &mined);
};


#endif

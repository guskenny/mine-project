/* Open Pit Mining(PCPSP)
 * Heuristic: construct an integer solution from a fractional solution
 * author: Davaa Baatar
 * date: 19.07.2016
*/
#ifndef __HEUR_FRAC_TO_INT_H__
#define __HEUR_FRAC_TO_INT_H__
#include "CumulativeModel.h"
//#include "typedef.h"
#include "BranchNode.h"
//#include "opt_distr.h" // what is this???
#include <algorithm>
#include "QOL/CpuTimer.h"
#include "compare.h"
#include "MineProblem.h"
// Construct an integer solution from fractional solution
// Ignore the branching history
class frac2Int{
public:
	frac2Int(CumulativeModel & prob, Sol_Real & sol, double time_limit = -1);
	~frac2Int();
	//enquiry
	bool getStatus();
	double getObjValue();
	double getCpuTime();
	const Sol_Int & getSolution();
	// method
	bool run();

private:
	// Input 
    	CumulativeModel & _prob;						// data
        
	const double _time_limit;						// time limit for the heuristic run; default -1 = no time limit;
	const double _disc_rate;						// discount factor
	
	// Output
	Sol_Int _sol;
	bool _status;									// return true if a feasible solution is constructed, 0 otherwise
	double _cpu_time;

	// internal
	double _current_factor;							
	CompareBlocks _compBlocks;						// 
	qol::CpuTimer _timer;							// timer
	std::vector<int> _Priority_list;       			// preference
	std::vector< std::vector<int> >	_start_list;	// period t -> index => block which starts at period t 
	void ignore(int it);
	void clean();
	int max_Priority();
};
#endif

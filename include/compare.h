/* Open Pit Mining(PCPSP)
 * Lexicographical order
 * author: Davaa Baatar
 * date: 19.07.2016
*/
#ifndef __COMPARE_H__
#define __COMPARE_H__
#include "CumulativeModel.h"
//#include "typedef.h"
#include <cfloat>
// Compare blocks lexicographically by :
	/*
	 *  - accumulated fractional values (up to the current period)
	 * 	- first period 
	 *  - Precedence relation
	 *  - maximum marginal profit(profit/resource usage rate), NOTE: for each block we might have different preference for destinations 
	 * 
	 * */
class CompareBlocks{ // use only for current period
public:
	CompareBlocks();
	CompareBlocks(CumulativeModel * prob, Sol_Real * sol);
	~CompareBlocks();
	
	void eval_nextPeriodPref(const std::vector<int> &A, const std::vector<int> &B);
	
	bool operator()(const int block_a, const int block_b);
	CompareBlocks & operator=(const CompareBlocks & rhs);
	void decrement_degree(int b);
	
private:
	// Input
	CumulativeModel * _prob; 				// problem data
	Sol_Real * _sol;						// fractional solution 
	
	// Internal
	int _period;							// current period (internal counter for the current active period)
	std::vector<double> _max_margProfit;	// maximum marginal profit
	std::vector<int> 	_start_period;      // starting period 
	std::vector<double> _preference;		// accumulated weights(update after each iteration)
	std::vector<int>    _access_degree;     // degree of accessibility
        std::vector<std::vector<int> >    _succs;   // successor blocsk for each block
};

class CompareDestinations{ 
public:
	CompareDestinations(const Block & b);
	~CompareDestinations();
	bool operator()(const int dest_a, const int dest_b);
private:

	// Input
	const Block & _block; 				// problem data
};

#endif


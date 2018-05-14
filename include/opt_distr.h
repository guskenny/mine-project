/* Open Pit Mining(PCPSP)
 * Exact algorithm: For a given period optimise the distribution of the allocated blocks 
 * author: Davaa Baatar
 * date: 19.07.2016
*/
#ifndef __OPT_DISTR_H__
#define __OPT_DISTR_H__
#include <util.h>
#include "CumulativeModel.h"
//#include "typedef.h"
#include "compare.h"
class OptDistr{
public:
	OptDistr(CumulativeModel & prob, Sol_Int & sol, const int t);
	~OptDistr();
	double getCpuTime();
	bool getStatus();
	bool solve(); // return 1 if solved, 0 otherwise
private:
	CumulativeModel & _prob;
	const int _period;
	Sol_Int & _sol;
	double _cpu_time;
	bool _status;
};
#endif

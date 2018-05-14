/* Open Pit Mining(PCPSP)
 * Lexicographical order
 * author: Davaa Baatar
 * date: 19.07.2016
*/
#include "compare.h"
CompareBlocks::CompareBlocks() : _prob(0), _sol(0), _period(0){}
CompareBlocks::CompareBlocks(CumulativeModel * prob, Sol_Real * sol) : _prob(prob), _sol(sol), _period(0) {
		
		const int nBlocks = (*_prob).getNBlock();
		const int t_max = (*_prob).getNPeriod();
		const int d_max = (*_prob).getnDestination();
		const int r_max = (*_prob).getnResources();
		// initialize the vectors 	_max_margProfit, _start_period, _preference
		_max_margProfit.resize(nBlocks, 0.0);
		_start_period.resize( nBlocks, t_max);
		_preference.resize(nBlocks, 0.0); // initialize
		_access_degree.resize(nBlocks, 0);
		_succs.resize(nBlocks);
		for(int b = 0; b< nBlocks; b++){
			const Block & block = (*_prob).getBlock(b);
			// get marginal rate
			double max_rate = DBL_MIN; 
			bool ghost = 1;
			for(int d=0; d<d_max; d++)
				for(int r=0; r<r_max; r++){
					if(block.getRCoef(d, r) > 1e-16){
						double rate = block.getProfit(d)/block.getRCoef(d, r);
						if(max_rate < rate) max_rate = rate;
						if(ghost) ghost = 0;
					}
			}
			if(ghost) max_rate = DBL_MAX;   //
			_max_margProfit[b] = max_rate; 

			// get start period
			for(int t=0; t<t_max; t++){
				if((*_sol).x[b][t] > 1e-15){
					_start_period[b] = t;
					break;
				}
			}
			// initialize the degree of accessibility
			_access_degree[b] = block.getNumPred();
			const std::vector<int> &preds=block.getPreds();
			for(auto p=preds.begin();p!=preds.end();++p)
			    _succs[*p].push_back(b);
		}
		
}
void CompareBlocks::decrement_degree(int b){ // be is removed
	for(auto s=_succs[b].begin();s!=_succs[b].end();++s)
	    -- _access_degree[*s];
}
CompareBlocks::~CompareBlocks(){}
void CompareBlocks::eval_nextPeriodPref(const std::vector<int> &A, const std::vector<int> &B){ 
		int n = (int) A.size();
		// increment weights for the remaing blocks in the heap
		for(int i=0; i<n; i++){
			int b = A[i];
			_preference[b] += _sol->x[b][_period]; 
		}
// increment weights for the entering blocks 
		n = (int) B.size();
		for(int i=0; i<n; i++){
			int b = B[i];
			_preference[b] += _sol->x[b][_period]; 
		}
		// increment period
		_period++;
}
bool CompareBlocks::operator()(const int a, const int b){ // block
    const int nBlocks = (*_prob).getNBlock();
    assert(0 <= a && 0 <= b && a < nBlocks && b < nBlocks);
    // check the degree of accessibility
    if(_access_degree[a] < _access_degree[b]) return 0;
    else if (_access_degree[a] > _access_degree[b]) return 1;
    else{	// check accumulated weights
	const double w_a = _preference[a];
	const double w_b = _preference[b];
	if(w_a > w_b + 1e-16) return 0;
	else if(w_a + 1e-16 < w_b) return 1;
	else{ //  w_a = w_b then check the earliest start
	    const int t_a = _start_period[a];
	    const int t_b = _start_period[b];
	    if(t_a < t_b) return 0; /* compare periods( t_a < t_b) then b is NOT a successor of a) */
	    else if (t_a > t_b) return 1;
	    else{ // starting at the same period check the precedence relation 
		const Block & block_a = (*_prob).getBlock(a);
		const Block & block_b = (*_prob).getBlock(b);
		if (block_a.isPred(b)) return 1; // true if b is pred of a
		else if (block_b.isPred(a)) return 0;
		else{ // NOT comparable check the marginal profit
		    double marg_a = _max_margProfit[a];
		    double marg_b = _max_margProfit[b];
		    if(marg_a > marg_b + eps_tol) return 0;
		}
	    }
	} 
    }
    return 1;
}
CompareBlocks & CompareBlocks::operator=(const CompareBlocks & rhs){
	_prob = rhs._prob; 		
	_sol = rhs._sol;		
	_period = rhs._period;
	_max_margProfit = rhs._max_margProfit;
	_start_period = rhs._start_period;      
	_preference = rhs._preference;
	_access_degree = rhs._access_degree;	
	_succs = rhs._succs;
	return *this;
}
CompareDestinations::CompareDestinations(const Block & b): _block(b){}
CompareDestinations::~CompareDestinations(){}
bool CompareDestinations::operator()(const int dest_a, const int dest_b){
	// order by profit( higher is better)
	double profit_a = _block.getProfit(dest_a);
	double profit_b = _block.getProfit(dest_b);
	if(profit_a < profit_b) return 1;
	else return 0;
}

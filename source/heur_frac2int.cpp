/* Open Pit Mining(PCPSP)
 * Heuristic: construct an integer solution from a fractional solution
 * author: Davaa Baatar
 * date: 19.07.2016
*/
#include "heur_frac2int.h"
#include "compare.h"
#include "opt_distr.h"

frac2Int::frac2Int(CumulativeModel & prob, Sol_Real & sol, double time_limit) :
    _prob(prob),  _time_limit(time_limit), _disc_rate(prob.getDiscountRate()),
    _sol(prob.getNBlock(),prob.getnDestination()),
    _status(0), _cpu_time(0.0), _compBlocks(&prob,&sol)  {
		
		const int nBlocks = _prob.getNBlock();
		const int t_max = _prob.getNPeriod();
		//const int d_max = _prob.getnDestination();
		//_sol = Sol_Int(nBlocks, d_max);		
		//_compBlocks = CompareBlocks(&prob, &sol);
		_current_factor = 1.0;
		_Priority_list.reserve(nBlocks); 
		_start_list = std::vector< std::vector<int> > (t_max, std::vector<int> ());
		for(int b = 0; b< nBlocks; b++){
			for(int t=0; t<t_max; t++){
				if( sol.x[b][t] > 1e-15){ 			// skip block if it is NOT in the fractional solution
					_start_list[t].push_back(b);	
					break;
				}
			}
		}
}
frac2Int::~frac2Int(){}
bool frac2Int::getStatus(){return _status;}
double frac2Int::getObjValue(){return _sol.obj;}
double frac2Int::getCpuTime(){return _cpu_time;}
const Sol_Int & frac2Int::getSolution(){ return _sol;}

void frac2Int::ignore(int it){ 
	_Priority_list[it] = -1; 
}
void frac2Int::clean(){
        //const int n = _Priority_list.size();
	std::vector<int>::iterator first =_Priority_list.begin();
	// push forward all non--negative indices
	for(std::vector<int>::iterator it = _Priority_list.begin(); it < _Priority_list.end(); it++){
		if(*it >-1){
			if(it != first) *first = *it;
			first++;
		}
	}
	_Priority_list.erase(first, _Priority_list.end());
}
int frac2Int::max_Priority(){ // return -1 if failed , otherwise index of the maximum element 
	int n = (int) _Priority_list.size();
	// get first active 
	int max_indx = -1;
	for(int i=0; i<n; i++){
		// check if it is active
		if(_Priority_list[i]>-1){
			max_indx = i;
			break;
		}
	}
	if(max_indx!=-1){
		for(int i = max_indx+1; i<n; i++){
			// check if it is active
			if(_Priority_list[i]>-1){
				if(_compBlocks(_Priority_list[max_indx], _Priority_list[i])) max_indx = i;
			}
		}	
	}	
	return max_indx;
}
bool frac2Int::run(){
    //const int nBlocks = _prob.getNBlock();
    const int t_max = _prob.getNPeriod();
    const int d_max = _prob.getnDestination();
    const int r_max = _prob.getnResources();
    bool stop = 0;
    int period = 0;		// internal counter for period
    _cpu_time = _timer.elapsedSeconds();  // start time
    double start_time =  _timer.elapsedSeconds();
    _sol.nT = 0;
    for(int b=0;b<_prob.getNBlock();++b) _sol.x[b] = t_max; // not mined
    while(!stop){
	// update accumulated weights
	_compBlocks.eval_nextPeriodPref(_Priority_list, _start_list[period]);
		
	// add new blocks 
	int n = (int) _start_list[period].size();
	for(int i=0; i<n; i++) 	_Priority_list.push_back(_start_list[period][i]); 
			
	int m = (int) _Priority_list.size();
	if(m>0){
					
	    // avialble resources for the period(Implemented for L type of constraints)
	    std::vector<double> resource_available(r_max, 0.0);
	    for(int r=0; r< r_max; r++) resource_available[r] = _prob.getLimit(r, period);
					
	    // add blocks one by one till first failure
	    for( int p=0; p<m; p++){
						
		/// get the max element	///
		int block_indx = max_Priority();
		if(block_indx == -1) continue; //skip 
		int block_id = _Priority_list[block_indx];
		const Block & block = _prob.getBlock(block_id);
		CompareDestinations _comp_d(block); // initialize compare
						
		// for each block order the destinations by the preference
		std::vector<int> _dest_pref;
		make_heap(_dest_pref.begin(), _dest_pref.end(), _comp_d); // empty
						
		for(int d=0; d<d_max; d++){
		    _dest_pref.push_back(d);
		    std::push_heap(_dest_pref.begin(), _dest_pref.end(), _comp_d);
		}
						
		// distribution 
		std::vector<int> _split(d_max, 0.0);
		double remaining = 1.0;
						
		for(int i =0; i<d_max; i++){
		    int d = _dest_pref.front();
		    std::pop_heap (_dest_pref.begin(),_dest_pref.end(), _comp_d); _dest_pref.pop_back();
		    double amount = 1.0;
		    for(int r=0; r< r_max; r++){
			double coef = block.getRCoef(d,r);
			if(coef>1e-11){
			    double limit = resource_available[r] / coef ;
			    if(amount > limit) amount = limit;
			}
		    } 
		    _split[d] = std::min(amount, remaining);
		    remaining -= _split[d];
		    if(remaining < 1e-8) break; 
		}
						
		if(remaining < 1e-8){// schedule the block
		    if(_sol.nT < period+1) _sol.nT = period + 1; // update number of active periods 
		    _sol.x[block_id] = period; 
		    for(int d=0; d<d_max; d++){
			// update obj value
			_sol.obj += _split[d]*block.getProfit(d)*_current_factor;
			// update available resources
			for(int r=0; r< r_max; r++) resource_available[r] -= _split[d]* block.getRCoef(d,r);
		    }
		    // mark the block for removal
		    ignore(block_indx); /// marked by -1
		    // update accessibilty degrees
		    _compBlocks.decrement_degree(block_id);
		}else{ // move to the next period
		    clean(); // clean the _Priority_list
		    break;
		}
	    }
	    clean();
	}
	// next period
	period++;
	_current_factor /= (1+_disc_rate);
			
	//check stopping criteria
	if(period == t_max) stop = true;
	if(_time_limit > 0  && _timer.elapsedSeconds()-start_time > _time_limit){
	    stop = true;
	}
    }
    // always have a sechedule, even if only some blocks included
    if(false && _Priority_list.size() > 0) _status = 0;  // failed to schedule
    else{	_status = 1;
	// reset the obj value
	_sol.obj = 0.0;
	for(int p = 0; p < _sol.nT; p++){ 
	    OptDistr reOpt(_prob, _sol, p);          
	    _status = reOpt.solve(); 	// re-optimize the schedule for each perido
	    if(!_status){
		std::cerr << "distribution optimisation failed for period " << p << std::endl;
		break; 	// return true if problem is solved
	    }
	}
    }
    _cpu_time = _timer.elapsedSeconds() - _cpu_time;
    // std::cout << "heur_frac2int: " << _status << ",  found " << getObjValue() << " in "
    // 	      << _cpu_time << " sec\n";
    return _status;
};
		

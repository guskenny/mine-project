#include "bz.h"
#include "heur_frac2int.h" // AE: need to re-add this

BZ::BZ( CumulativeModel &prob, MaxClosureFactory mcfactory, double timeLimit,const BranchNode_info *branch): 
  _nIter(0), _prob(prob), 
  _branch(prob.getNBlock(),prob.getNPeriod(),prob.getnDestination(),prob.getnResources()),
  _status(0), _timelimit(timeLimit), 
  _ub(0.0), _obj_const(0.0), _obj_val(0.0), _cpu_mc(0.0), _cpu_rlp(0.0),
  _Cpart(prob,branch) {
  
    if( branch){
	_branch = *branch; // copy branch information
	std::cout << "Branch has " << _branch.fixedCnt() << " fixed\n";
    }

	const int nBlocks = _prob.getNBlock();
	const int t_max = _prob.getNPeriod();
	const int d_max = _prob.getnDestination(); 
	const int r_max = _prob.getnResources();
	
	int n = t_max*r_max;
	_lambda = std::vector<double> (n, 0.0);
	
	n = nBlocks*t_max*d_max;
	_supply = std::vector<double> (n, 0.0);
	_nVertices = nBlocks*t_max*d_max ;
	_maxCut = std::vector<int> (_nVertices, 0);
	_sol_z = std::vector<double> (n, 0.0); 
	// set up the time limit
	if(_timelimit <= 0) {
		_timer.setTimeLimit(1e20);
    }else _timer.setTimeLimit(_timelimit);
    
    _cpu_network = _timer.elapsedSeconds();
    std::map<int,int> fixed;
    _prob.BN_info_to_map(_branch,fixed);
    std::cout << "Num variables fixed: " << fixed.size() << std::endl;
    _mc = mcfactory(_prob.graph,fixed); 
	
	_cpu_network = _timer.elapsedSeconds() - _cpu_network;
	std::cout<<(boost::format("Network is constructed in %7.2f sec.") %_cpu_network)<<std::endl;
	_cpu_time = _timer.elapsedSeconds();
	
}
BZ::~BZ(){
	if(_mc!=NULL) delete _mc;
	_mc = NULL;
}
bool BZ::solve(MineProblem *best){
	
	std::cout << "Branch has " << _branch.fixedCnt() << " fixed\n";

	std::cout << "----------------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "Bienstock-Zuckerberg Algorithm (Initial Version) " << std::endl;
	std::cout << "----------------------------------------------------------------------------------------------------------" << std::endl;
	std::cout<<(boost::format("%10s   %12s   %9s   %12s   %9s   %9s   %11s   %12s") % " Iteration" % "LR" %"cpu" % "RLP" %"cpu" %"Gap %" %"CPU"%"Heur")<<std::endl;
	std::cout << "----------------------------------------------------------------------------------------------------------" << std::endl;
	
	// start solving
	double tmp, tmp1;
	
	/// Bienstock-Zuckerberg Algorithm
	while(!_status && ! _timer.timeLimitReached()){
		_nIter++;
		tmp = _timer.elapsedSeconds(); 
		// evaluate the constant term and coefficients of the corresp. obj. function
		_obj_const = _prob.calcRedCost(_lambda, _supply);
		_ub = _obj_const;
		std::cout<<(boost::format("%5s%5d") %"Iter." %_nIter);
		
		(*_mc).setProfit(_supply);
		_status = (*_mc).solve(); // return non-zero if not solved
		tmp1 = _timer.elapsedSeconds() - tmp;
		_cpu_mc += tmp1;
		if(!_status){// solved
			_status = (*_mc).getClosure(_maxCut); // return zero if failed
			if(_status){ // obtained maxcut
				_ub += (*_mc).calcProfit(_supply, _maxCut);
				// update partitioning
				_status = updatePartitioning(); // return 1 if stop criteria holds
				std::cout << (boost::format("   %12.0f   %9.2f") %_ub %tmp1);
				if(_status) std::cout<<"\t Stop, optimal "<<std::endl; 
			}else{
				std::cerr<<"Error obtaining network solution (node values) ... " << _status << std::endl;
				break;
			}
		}else{
			 std::cerr<<"Abort : Could NOT solve the MaxClosure Problem! "<<std::endl;
			 break;
		} 
		
		if(!_status){// solve restricted problem with new additional constraints defined w.r.t. Cpart
			// convert _lambda to 2D vector
			const int t_max = _prob.getNPeriod();
			const int r_max = _prob.getnResources();
			std::vector<std::vector<double> > lambda(t_max, std::vector<double> (r_max, 0.0));
			tmp1 = _timer.elapsedSeconds();
			ResLP lp_res(_prob, _Cpart);  // each time new lp is constructed
			/* optimal solution, lambda and stopping criteria are updated automatically. Ugly but less memory usage ! */
			lp_res.solve(_sol_z,lambda);
			tmp1 = _timer.elapsedSeconds() - tmp1;
			tmp = _timer.elapsedSeconds() - tmp;
			_cpu_rlp += tmp1;
			if(lp_res.getStatus()){ 
				_obj_val = lp_res.getObjVal();
				restore_Multipliers(lambda); //	lp_res.solve_LargeLP(); 
				std::cout << (boost::format("   %12.0f   %9.2f   %9.6f   %11.2f")
					      %_obj_val %tmp1 %((_ub - _obj_val)/_obj_val) %tmp
				    );
			}else{
				std::cout<<"   Abort: Restricted LP not solved!\n";
				break;
			}
			if(best){			// run heuristic
				Sol_Real fracSol(_prob.getNBlock(),t_max,_prob.getnDestination());
				// get fractional solution length nB*nT*nD
				_prob.storeSoln(fracSol,_sol_z);
				frac2Int heur(_prob,fracSol); // don't worry about time limit
				heur.run();
				if(heur.getStatus()){
				    best->update_sol(heur.getSolution());
				    std::cout << boost::format("   %12.0f")%heur.getObjValue()
 					      << (best->get_obj()==heur.getObjValue() ? '*' : ' ');
				    
				}else
				    std::cout << "failed";
			}
			std::cout << std::endl;
		}
	}
	if(best){
	    Sol_Real sol(_prob.getNBlock(),_prob.getNPeriod(),_prob.getnDestination());
	    // get fractional solution length nB*nT*nD
	    _prob.storeSoln(sol,getSolution());
	    sol.obj=_obj_val;
	    best->update_sol(sol);
	}
	_cpu_time = _timer.elapsedSeconds() - _cpu_time;
	std::cout << "----------------------------------------------------------------------------------------------------------" << std::endl;
	return _status; 
}

void BZ::reindex_Multipliers(std::vector<std::vector<double> > & lambda){
	const int t_max = _prob.getNPeriod();
	const int r_max = _prob.getnResources();
	int index = 0;
	for(int t=0; t<t_max; t++)
		for(int r=0; r<r_max; r++){
				lambda[t][r] = _lambda[index];
				index++;
	}
}

void BZ::restore_Multipliers(const std::vector<std::vector<double> > & lambda){
	const int t_max = _prob.getNPeriod();
	const int r_max = _prob.getnResources();
	int index = 0;
	for(int t=0; t<t_max; t++)
		for(int r=0; r<r_max; r++){
				_lambda[index] = lambda[t][r];
				index++;
	}
}
int BZ::updatePartitioning(){
    
    return _Cpart.split(_maxCut)==0;
}
double BZ::getCpuTime(){return _cpu_time;}
double BZ::getCpuTime_MC(){return _cpu_mc;}
double BZ::getCpuTime_RLP(){return _cpu_rlp;}
double BZ::getCpuTime_Network(){return _cpu_network;}
int BZ::get_nIter(){return _nIter;}
int BZ::getStatus(){return _status;}
double BZ::getObjValue(){return _obj_val;}
double BZ::getUB(){return _ub;}
std::vector<double> & BZ::getSolution(){return _sol_z;}

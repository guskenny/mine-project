/* Open Pit Mining(PCPSP)
 * Exact algorithm: For a given period optimise the distribution of the allocated blocks 
 * author: Davaa Baatar
 * date: 19.07.2016
*/
#include <ilcplex/ilocplex.h>
#include "opt_distr.h"
#include <math.h> 
ILOSTLBEGIN
typedef IloArray< IloNumVarArray > IloNumVarArray2;

OptDistr::OptDistr(CumulativeModel & prob, Sol_Int & sol, const int t) : _prob(prob), _sol(sol), _period(t), _cpu_time(0.0), _status(0) {

// print input solution
/*	const int nBlocks =_prob.getNBlock();
	const int d_max =_prob.getnDestination();
	std::cout<<"period:\t"<<	_period<<std::endl;
	for(int b=0; b<nBlocks; b++){
		std::cout<<(boost::format("%6i") %b);
		if(sol.X[b]> -1e-10){
				std::cout<<(boost::format("  %8i  %2i")%1 %sol.X[b]);
					for(int d = 0; d<d_max; d++){
							std::cout<<(boost::format("  %8.6f")%sol.Y[b][d]);
					}
					std::cout<<std::endl;
				}else{ 
					std::cout<<(boost::format("  %8i")%0)<< std::endl;
				}
	}
	*/
}
OptDistr::~OptDistr(){}
double OptDistr::getCpuTime() {return _cpu_time;}
bool OptDistr::getStatus(){return _status;}
bool OptDistr::solve(){ // return 1 if solved, 0 otherwise
IloEnv env;
	IloTimer timer(env);
	timer.start();
	try{
		const int r_max =_prob.getnResources();
		const int nBlocks =_prob.getNBlock();
		//const int t_max =_prob.getNPeriod();
		const int d_max =_prob.getnDestination();
		std::vector<Block> * blocks = _prob.getBlock();
	// mappings
		std::vector<int> active_blocks;
		active_blocks.reserve(nBlocks);
		for(int b=0; b<nBlocks; b++){
			if(_sol.x[b] == _period) active_blocks.push_back(b); 
		}
	
	
	// formulate the LP in terms of indices of the active blocks in the vector
		const int nBlocks_active = (int) active_blocks.size();
		//const int nVar = d_max*nBlocks_active;
	
	// decision variables
		IloNumVarArray2 Yvar(env);  // block -> destination => value(%)
		for(int b=0; b<nBlocks_active; b++)
			 Yvar.add(IloNumVarArray(env, d_max, 0, 1, ILOFLOAT));
		 
		 // type of constraints
		 
		 /* Distribution */
		 IloRangeArray SumDest(env);
		 for(int b=0; b<nBlocks_active; b++){
			IloExpr expr(env);
			for(int d=0; d<d_max; d++)
					expr += Yvar[b][d];
			SumDest.add(IloRange(env, 1, expr, 1));
		}
		 
		 /* Resource constraints (10) */
		 IloRangeArray Resource(env);
		 for(int r=0; r<r_max; r++){
				char cType = _prob.getResConstrType(r, _period);
				IloExpr expr(env);
				for(int b=0; b<nBlocks_active; b++){
					int id = active_blocks[b];  // mapping back to actual block
					for(int d=0; d<d_max; d++){
						double coef = (*blocks)[id].getRCoef(d,r);
						expr += coef*Yvar[b][d];
					}
				}
				if(cType == 'L'){ 
					Resource.add(IloRange(env, -IloInfinity, expr, _prob.getLimit(r, _period)));
				}else if(cType == 'R'){
						Resource.add(IloRange(env, _prob.getLimit(r, _period), expr, IloInfinity));
					 }else if(cType == 'I'){
						 Resource.add(IloRange(env, _prob.getLimit(r, _period), expr, _prob.getLimit(r, _period, 1)));
					 }
		 }
		 
		 // objective function
		 IloExpr exprObj(env);
		 IloNum discount_factor = 1/IloPower(1+_prob.getDiscountRate(), _period);
		 for(int b=0; b<nBlocks_active; b++){
			int id = active_blocks[b];  // mapping back
			for(int d=0; d<d_max; d++){
				double coef = (*blocks)[id].getProfit(d)*discount_factor;
				exprObj += coef*Yvar[b][d];
			}
		 }
		 IloObjective Obj=IloMaximize(env, exprObj);
		 // cplex model 
		 IloModel distr_opt(env);
		 distr_opt.add(SumDest);
		 distr_opt.add(Resource);
		 distr_opt.add(Obj);
		 IloCplex distr_opt_cpx(distr_opt);
//		 distr_opt_cpx.exportModel("tyni.lp");
		 distr_opt_cpx.setParam(IloCplex::ClockType,1);
		 distr_opt_cpx.setOut(env.getNullStream());
		 // solve & retrieve related info
		 timer.start();
		 distr_opt_cpx.solve();
		 timer.stop();
		 
		 // extruct solution	
		 for(int b=0; b<nBlocks_active; b++){
			 IloNumArray sol(env);
			 distr_opt_cpx.getValues(sol, Yvar[b]);	
			 int id = active_blocks[b];  // mapping back
			 for(int d=0; d<d_max; d++){
				_sol.y[id][d] = sol[d];
				// update obj	
				double coef = (*blocks)[id].getProfit(d)*discount_factor;
				_sol.obj += coef*_sol.y[id][d]; 
			}
		 }
// print solution on screen
/*		 std::cout<<"\t New Solution"<<std::endl;
			 
		for(int b=0; b<nBlocks; b++){
			std::cout<<(boost::format("%6i") %b);
			if(_sol.X[b]> -1e-10){
					std::cout<<(boost::format("  %8i  %2i")%1 %_sol.X[b]);
						for(int d = 0; d<d_max; d++){
								std::cout<<(boost::format("  %8.6f")%_sol.Y[b][d]);
						}
						std::cout<<std::endl;
					}else{ 
						std::cout<<(boost::format("  %8i")%0)<< std::endl;
					}
		}*/
		 _status = 1;
		// free memory
		for(int b=0; b<nBlocks_active; b++){
				Yvar[b].endElements();
				Yvar[b].end();
		}
		Yvar.end(); 

	}
	catch (IloException & ex) { cerr << "Error: " << ex << endl; }
	catch (...) { cerr << "Unknown Error Occured" << endl;}
	timer.stop();
	_cpu_time = timer.getTime();
	env.end();
	return _status;
}	
	

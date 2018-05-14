/* LP Relaxation: Open Pit Mining(PCPSP)
 * author: Davaa Baatar
 * date: 19.06.2016
*/
#include <rlp.h>
#include <ilcplex/ilocplex.h>
#include "QOL/QolMIP.h"
#include "QOL/QolColFormulation.h"
#include "QOL/CplexFormulation.h"
#include "QOL/GurobiFormulation.h"
#include "QOL/CpuTimer.h"
#include "QOL/QolParams.h"
#include "QOL/QolUtil.h"


ILOSTLBEGIN
typedef IloArray< IloRangeArray > IloRangeArray2;

void Partition::setFixed(const Daten &model,const BranchNode_info &branch,bool splitAll)
{
    const int nBlocks =model.getNBlock();
    const int t_max =model.getNPeriod();
    const int d_max =model.getnDestination();
    
    resize(3);
    for(auto p=begin();p!=end();++p) p->clear(); // clear each partition
    (*this)[0].reserve(nBlocks*t_max*d_max);
    fixed.resize(3);
    fixed[0]=-1;
    fixed[1]=0;
    fixed[2]=1;
    for(int a=0; a<nBlocks; a++){
	int t;
	for(t=0;t<t_max;++t)
	    for(int d=0;d<d_max-1;++d)
		(*this)[0].push_back(convert_triplet_to_index(model, a, t, d));
	for(t=0; t<branch.time[a][0]; t++){
	    int index_a = convert_triplet_to_index(model, a, t, d_max-1);
	    (*this)[1].push_back(index_a);
	}
	for(; t<branch.time[a][1]; t++){
	    int index_a = convert_triplet_to_index(model, a, t, d_max-1);
	    (*this)[0].push_back(index_a);
	}
	// FIXME: should check branch.mine_block
	for(; t<t_max; t++){
	    int index_a = convert_triplet_to_index(model, a, t, d_max-1);
	    (*this)[2].push_back(index_a);
	}
    }
    for(int i=2;i>=0;--i)
	if( (*this)[i].empty() ){
	    erase(begin()+i);
	    fixed.erase(fixed.begin()+i);
	}
    if(splitAll && fixed[0]==-1 && begin()->size() > 1){
	size_t end=size()-1;
	resize(size()+begin()->size()-1);
	for(size_t i=1;i<begin()->size();++i)
	    (begin()+end+i)->push_back( (*begin())[i]);
	begin()->resize(1);
	fixed.resize(size(),-1);
    }
    
}


int Partition::split(const std::vector<int> &cut)
{
    size_t len= size();
    for(size_t i=0; i< len; i++){
	if( isFixed(i) ){
	    for(int v : (*this)[i]){
		if( cut[v] != fixed[i] )
		    std::cerr <<"ERROR cut["<<v
			      <<"] should be " << fixed[i]
			      << " not " << cut[v] << std::endl;
	    }
	    continue; // cut should never split this set
	}
	const std::vector<int> &part=(*this)[i];
	int m = (int) part.size();
	if(m==1) continue;
	std::vector<int> I; I.reserve(m);
	std::vector<int> O; O.reserve(m);
	for(int j: part){
	    if(cut[j]>=1) I.push_back(j);
	    else O.push_back(j);
	}
	if(I.size()>0 && O.size() > 0){
	    (*this)[i]=I;
	    push_back(O);
	}
    }
    if(size() > len) fixed.resize(size(),-1); // fill with non-fixed
    return int(size()-len);
}

// Simple implementation: just split based on integer values
// Should be a bit more sophisticated to enable fractional y's to be
// reproduced correctly
int Partition::split(const CumulativeModel &model,const Sol_Int &soln)
{
    std::vector<int> y(n,0.0);
    model.getSoln(soln,y);    
    std::vector<int> s(n,0);
    for(int i=0;i<n;++i)
	s[i] = int(y[i] > 0.5);	// round up
    return split(s);
}



ResLP::ResLP(CumulativeModel & data, const Partition& Cpart):
	_data(data), _Cpart(Cpart){}
ResLP::~ResLP(){}
double ResLP::getObjVal(){return _objVal;}
double ResLP::getTime(){return _t_cpu;}

void ResLP::solve(std::vector<double> & _sol,std::vector<std::vector<double> > & _pi){
    IloEnv env;
    IloTimer timer(env);
    timer.start();
    try{
	const int r_max =_data.getnResources();
	const int nBlocks =_data.getNBlock();
	const int t_max =_data.getNPeriod();
	const int d_max =_data.getnDestination();
	const int nVar = (int) _Cpart.size();// + dummy fixed to zero/one
	const int nPart = (int) _Cpart.size();
	const int nZ = nBlocks*t_max*d_max;
	std::vector<int> part_mapping(nZ, -1);
	IloNumVarArray Z(env, nVar, 0, 1, ILOFLOAT);
	for(int v=0;v<nPart;++v){
	    if(_Cpart.isFixed(v)){
		Z[v].setUB(_Cpart.fixedVal(v));
		Z[v].setLB(_Cpart.fixedVal(v));
	    }
	    for(int index : _Cpart[v])
		part_mapping[index] = v;
	}
	std::vector<Block> * blocks =_data.getBlock();
	/* CPLEX model */
		
	IloRangeArray Precedence(env); /// prec constraints
	{   // (12)
	    std::vector< std::vector<int> > ind(nVar, std::vector<int>(nVar,0));
	    for(int a=0; a<nBlocks; a++)
		for(int t=0; t<t_max; t++){
		    std::vector<int> * pred = (*blocks)[a].getPreds();
		    for(int i=0; i< (*blocks)[a].getNumPred(); i++){
			int b = (*pred)[i];
			int index_a = convert_triplet_to_index(_data, a, t, d_max-1);
			int index_b = convert_triplet_to_index(_data, b, t, d_max-1);
			index_a = part_mapping[index_a]; // successor
			index_b = part_mapping[index_b]; // predecessor
			if(index_a != index_b){
			    if(ind[index_a][index_b]==0){
				ind[index_a][index_b]=1;
				IloExpr expr(env);
				expr = Z[index_a] - Z[index_b];
				Precedence.add(IloRange(env, -1, expr, 0));
			    }
			}
		    }
		}
		
	    // (14) 
	    for(int b=0; b<nBlocks; b++)
		for(int t=0; t<t_max; t++)
		    for(int d=1; d<d_max; d++){
			int index_l = convert_triplet_to_index(_data, b, t, d-1);
			int index_r = convert_triplet_to_index(_data, b, t, d);
			index_l = part_mapping[index_l];
			index_r = part_mapping[index_r];
						
			if(index_l != index_r){
			    if(ind[index_l][index_r]==0){
				ind[index_l][index_r]=1;
				IloExpr expr(env);
				expr = Z[index_l] - Z[index_r];
				Precedence.add(IloRange(env, -1, expr, 0));
			    }
			}
		    }
	    // (15)
	    for(int b=0; b<nBlocks; b++)
		for(int t=1; t<t_max; t++){
		    int index_l = convert_triplet_to_index(_data, b, t-1, d_max-1);
		    int index_r = convert_triplet_to_index(_data, b, t, 0);
		    index_l = part_mapping[index_l];
		    index_r = part_mapping[index_r];
					
		    if(index_l != index_r){
			if(ind[index_l][index_r]==0){
			    ind[index_l][index_r]=1;
			    IloExpr expr(env);
			    expr = Z[index_l] - Z[index_r];
			    Precedence.add(IloRange(env, -1, expr, 0));
			}
		    }
		}
	}		
	    IloRangeArray2 Resource(env);
	    // (16) add constraints in the same order as lambda is constructed
	    for(int t=0; t<t_max; t++){
		Resource.add(IloRangeArray(env));
		for( int r=0; r<r_max; r++){
		    std::vector<double> coef(nVar, 0.0);
		    std::vector<int> ind(nVar, 0);/// could be used limits, instead
		    for(int b=0; b<nBlocks; b++){
			int index;
			/// q_{1,r}z_{bt1}
			index = convert_triplet_to_index(_data, b, t, 0);
			index = part_mapping[index];
			if(ind[index]==0) ind[index] = 1; 
			coef[index] += (*blocks)[b].getRCoef(0,r);
						
			/// sum_{d=1}^{d_max} q_{dr}*(z_{btd}-z_{bt,d-1}) 
			for(int d=1; d<d_max; d++){
			    /// z_{btd}
			    index = convert_triplet_to_index(_data, b, t, d);
			    index = part_mapping[index];
			    if(ind[index]==0) ind[index] = 1; 
			    coef[index] += (*blocks)[b].getRCoef(d,r);
			    /// -z_{bt,d-1}
			    index = convert_triplet_to_index(_data, b, t, d-1);
			    index = part_mapping[index];
			    if(ind[index]==0) ind[index] = 1; 
			    coef[index] -= (*blocks)[b].getRCoef(d,r);
			}
			if(t>0){
			    /// -q_{1r}*z_{b,t-1,d_max}
			    index = convert_triplet_to_index(_data, b, t-1, d_max-1);
			    index = part_mapping[index];
			    if(ind[index]==0) ind[index] = 1; 
			    coef[index] -= (*blocks)[b].getRCoef(0,r);
			}
		    }
				
		    IloExpr expr(env);
		    for(int i=0; i<nVar; i++)
			if(ind[i]) expr += coef[i]*Z[i];
		    Resource[t].add(IloRange(env, -IloInfinity, expr,_data.getLimit(r, t)));  ///only type L
		}
	    }
				
	    IloExpr exprObj(env);
	    {
		std::vector<double> coef(nVar, 0.0);
		const double rate =_data.getDiscountRate();
		for(int b=0; b<nBlocks; b++){
		    double cur_rate = 1.0;
		    for(int t=0; t<t_max; t++){
			for(int d=0; d<d_max; d++){
			    int index = convert_triplet_to_index(_data, b, t, d);
			    index = part_mapping[index];
			    coef[index] += (*blocks)[b].getProfit(d)*cur_rate; 
			    if(d<d_max-1) coef[index] -= (*blocks)[b].getProfit(d+1)*cur_rate; 
			    else
				if(t<t_max-1) coef[index] -= (*blocks)[b].getProfit(0)*cur_rate/(1+rate); 
			}
			cur_rate /= (1+rate);
		    }
		}
			
		for(int i=0; i<nVar; i++) exprObj += coef[i]*Z[i];
	    }
	    IloObjective Obj=IloMaximize(env, exprObj);
	    // cplex model 
	    IloModel restLP(env);
	    restLP.add(Precedence);
	    for(int t=0; t<t_max; t++)
		restLP.add(Resource[t]);
	    restLP.add(Obj);
	    //restLP.add( Z[ZeroVar] <= Z[OneVar]); // make sure variables are added
		 
	    IloCplex rlp_cpx(restLP);
	    rlp_cpx.setOut(env.getNullStream());
	    rlp_cpx.setParam(IloCplex::ClockType,1);
	    rlp_cpx.solve();
	    if(rlp_cpx.getStatus() == IloAlgorithm::Optimal) _status = 1;
	    else _status = 0; 
	    if(_status){
		IloNumArray sol_rlp(env);
		rlp_cpx.getValues(sol_rlp, Z);
		for(int i=0; i<nZ; i++){
		    _sol[i] = sol_rlp[ part_mapping[i] ];
		}
		_objVal = (double) rlp_cpx.getObjValue();
		for(int t=0; t<t_max; t++){
		    rlp_cpx.getDuals(sol_rlp, Resource[t]);
		    for( int r=0; r<r_max; r++)	_pi[t][r]= sol_rlp[r];
		}
	    }
	    // free memory
	    for(int t=0; t<t_max; t++){
		Resource[t].endElements();
		Resource[t].end();
	    }
	    Resource.end(); 
    }
    catch (IloException & ex) { cerr << "Error: " << ex << endl; }
    catch (...) { cerr << "Unknown Error Occured" << endl;}
    timer.stop();
    _t_cpu = timer.getTime();
    env.end();
} // end solve()


void ResLP::solveExact(MineProblem &prob, Sol_Int &best,double timeLimit)
// This implementation is based on solve() but
// (a) Uses QOL to make switching solvers easier (defaults to GUROBI)
// (b) Enforces block variables to be binary
// (c) Extracts an integer solution and does not produce any dual information    
{
    
    const qol::FormulationType solverType=qol::GUROBI; 
    qol::MIPSolver *mip =0;
    qol::Parameters param;
    qol::CpuTimer timer;
    param.setParamVal(qol::VERBOSITY,1);
    param.setParamVal(qol::WALLTIME,timeLimit);
    qol::Status status = qol::ABORTED;
    switch(solverType){
	case qol::GUROBI: mip = new qol::GurobiFormulation(); break;
	case qol::CPLEX: mip = new qol::CplexFormulation(); break;
	default:
	    std::cerr << "QOL Solver type not supported!";
    }
    try{
	const int r_max =_data.getnResources();
	const int nBlocks =_data.getNBlock();
	const int t_max =_data.getNPeriod();
	const int d_max =_data.getnDestination();
	const int nVar = (int) _Cpart.size();// + dummy fixed to zero/one
	const int nZ = nBlocks*t_max*d_max;
	std::vector<int> part_mapping(nZ, -1);
	std::vector<qol::Variable> Z=mip->addVarVec(nVar);
	std::vector<char> binary(nZ,qol::Variable::CONTINUOUS);
	for(int v=0;v<nVar;++v){
	    mip->setVarName(Z[v],boost::str(boost::format("z%03d")%v));
	    //char binary=qol::Variable::CONTINUOUS;	// does part include at least one binary
	    for(int index : _Cpart[v]){
		part_mapping[index] = v;
		int b,t,d;
		convert_index_to_triplet(_data,v,b,t,d);
		if(d == d_max-1) binary[v] = qol::Variable::BINARY;
	    }
	    if(_Cpart.isFixed(v)){
		mip->setVarLB(Z[v],_Cpart.fixedVal(v));
		mip->setVarUB(Z[v],_Cpart.fixedVal(v));
	    }else{
		mip->setVarLB(Z[v],0);
		mip->setVarUB(Z[v],1);
		mip->setVarType(Z[v],binary[v]); // if it is binary
	    }
	}
	std::vector<Block> * blocks =_data.getBlock();
	/* CPLEX model */
	
	//IloRangeArray Precedence(env); /// prec constraints
	{   // (12)
	    std::vector< std::vector<int> > ind(nVar, std::vector<int>(nVar,0));
	    for(int a=0; a<nBlocks; a++){
		std::vector<int> * pred = (*blocks)[a].getPreds();
		for(int t=0; t<t_max; t++)
		    for(int b : (*pred) ){
			int index_a = convert_triplet_to_index(_data, a, t, d_max-1);
			int index_b = convert_triplet_to_index(_data, b, t, d_max-1);
			index_a = part_mapping[index_a]; // successor
			index_b = part_mapping[index_b]; // predecessor
			if(index_a != index_b && ind[index_a][index_b]==0){
			    ind[index_a][index_b]=1;
			    mip->addConstraint(Z[index_a] <= Z[index_b]);
			    //IloExpr expr(env); expr = Z[index_a] - Z[index_b];
			    //Precedence.add(IloRange(env, -1, expr, 0));
			}
		    }
	    }
	    
	    // (14) 
	    for(int b=0; b<nBlocks; b++)
		for(int t=0; t<t_max; t++)
		    for(int d=1; d<d_max; d++){
			int index_l = convert_triplet_to_index(_data, b, t, d-1);
			int index_r = convert_triplet_to_index(_data, b, t, d);
			index_l = part_mapping[index_l];
			index_r = part_mapping[index_r];						
			if(index_l != index_r && ind[index_l][index_r]==0){
			    ind[index_l][index_r]=1;
			    mip->addConstraint(Z[index_l] <= Z[index_r]);
			}
		    }
	    // (15)
	    for(int b=0; b<nBlocks; b++)
		for(int t=1; t<t_max; t++){
		    int index_l = convert_triplet_to_index(_data, b, t-1, d_max-1);
		    int index_r = convert_triplet_to_index(_data, b, t, 0);
		    index_l = part_mapping[index_l];
		    index_r = part_mapping[index_r];
		    
		    if(index_l != index_r && ind[index_l][index_r]==0){
			ind[index_l][index_r]=1;
			mip->addConstraint(Z[index_l] <= Z[index_r]);
		    }
		}
	} // end precedence constraints
	    
	std::vector<std::vector<qol::Constraint> >  Resource(t_max);
	// (16) add constraints in the same order as lambda is constructed
	for(int t=0; t<t_max; t++){
	    for( int r=0; r<r_max; r++){
		std::vector<double> coef(nVar, 0.0);
		std::vector<int> ind(nVar, 0);/// could be used limits, instead
		for(int b=0; b<nBlocks; b++){
		    int index;
		    /// q_{1,r}z_{bt1}
		    index = convert_triplet_to_index(_data, b, t, 0);
		    index = part_mapping[index];
		    if(ind[index]==0) ind[index] = 1; 
		    coef[index] += (*blocks)[b].getRCoef(0,r);
		    
		    /// sum_{d=1}^{d_max} q_{dr}*(z_{btd}-z_{bt,d-1}) 
		    for(int d=1; d<d_max; d++){
			/// z_{btd}
			index = convert_triplet_to_index(_data, b, t, d);
			index = part_mapping[index];
			if(ind[index]==0) ind[index] = 1; 
			coef[index] += (*blocks)[b].getRCoef(d,r);
			/// -z_{bt,d-1}
			index = convert_triplet_to_index(_data, b, t, d-1);
			index = part_mapping[index];
			if(ind[index]==0) ind[index] = 1; 
			coef[index] -= (*blocks)[b].getRCoef(d,r);
		    }
		    if(t>0){
			/// -q_{1r}*z_{b,t-1,d_max}
			index = convert_triplet_to_index(_data, b, t-1, d_max-1);
			index = part_mapping[index];
			if(ind[index]==0) ind[index] = 1; 
			coef[index] -= (*blocks)[b].getRCoef(0,r);
		    }
		}
		
		qol::Expression expr;
		for(int i=0; i<nVar; i++)
		    if(ind[i]) expr += coef[i]*Z[i];
		Resource[t].push_back(mip->addConstraint(expr <= _data.getLimit(r, t)));  ///only type L
	    }
	} // endloop over t
				
	{
	    std::vector<double> coef(nVar, 0.0);
	    const double rate =_data.getDiscountRate();
	    for(int b=0; b<nBlocks; b++){
		    double cur_rate = 1.0;
		    for(int t=0; t<t_max; t++){
			for(int d=0; d<d_max; d++){
			    int index = convert_triplet_to_index(_data, b, t, d);
			    index = part_mapping[index];
			    coef[index] += (*blocks)[b].getProfit(d)*cur_rate; 
			    if(d<d_max-1) coef[index] -= (*blocks)[b].getProfit(d+1)*cur_rate; 
			    else
				if(t<t_max-1) coef[index] -= (*blocks)[b].getProfit(0)*cur_rate/(1+rate); 
			}
			cur_rate /= (1+rate);
		    }
		}
			
	    for(int i=0; i<nVar; i++) mip->setObjCoeff(Z[i],-coef[i]); // minimise -ve
	}
	mip->setParameters(param);
	//TODO: Input current best solution as heuristic bound for solve
	//      Add a MIP call-back to communicate improved solutions whenever
	//      they are found.
	status = mip->optimize();
	if(status== qol::HEURISTIC || status==qol::OPTIMAL) { // extract solution
	    _status = 1;
	    std::vector<double> sln(_data.graph.getNumNodes(),0.0);
	    size_t cnt=0;
	    double profit=0.0;
	    for(int v=0;v<nVar;++v){
		if(_Cpart.isFixed(v)){
		    if( _Cpart.fixedVal(v) == 1) cnt += _Cpart[v].size();
		    for(int index : _Cpart[v]){
			sln[index]=_Cpart.fixedVal(v);
			profit += sln[index]*_data.getProfit(index);
		    }
		}else{
		    double zval = mip->getPrimal(Z[v]);
		    if(binary[v]!=qol::Variable::CONTINUOUS)
			zval = floor(zval+0.5); // make sure we are integer
		    if( zval > 0.5) cnt +=  _Cpart[v].size();
		    for(int index : _Cpart[v]){
			sln[index] = zval;
			profit += zval*_data.getProfit(index);
		    }
		}
	    }
	    _data.storeSoln(best,sln);
	    prob.update_sol(best);
	    std::cout << "\tResLP completed with solution " << best.obj
		      << " = " << profit
		      << " using " << best.nT << " periods" << std::endl;
	}else{
	    std::cerr << "\tERROR: cannot solve, status = " << status << std::endl;
	    mip->writeLP("infeasible.lp");
	    std::cerr << "\t Wrote infeasible.lp\n";
	}
    }
    catch (qol::Exception & ex) { cerr << "Error: " << ex.what() << std::endl; }
    catch (...) { cerr << "Unknown Error Occured" << endl;}
    _t_cpu = timer.elapsedSeconds(); // CPU time
} // end solveExact()


void ResLP::solve_LargeLP(std::vector<double> & _sol,std::vector<std::vector<double> > & _pi){
	IloEnv env;
	IloTimer timer(env);
	timer.start();
	try{
		const int r_max =_data.getnResources();
		const int nBlocks =_data.getNBlock();
		const int t_max =_data.getNPeriod();
		const int d_max =_data.getnDestination();
		const int nZ = nBlocks*t_max*d_max;
				
		std::vector<Block> * blocks =_data.getBlock();
		
		/* CPLEX model */
		
		IloNumVarArray Z(env, nZ, 0, 1, ILOFLOAT);
		
		IloRangeArray Precedence(env); /// prec constraints
		   // (12)
		for(int a=0; a<nBlocks; a++)
			for(int t=0; t<t_max; t++){
				std::vector<int> * pred = (*blocks)[a].getPreds();
				for(int i=0; i< (*blocks)[a].getNumPred(); i++){
					int b = (*pred)[i];
					int index_a = convert_triplet_to_index(_data, a, t, d_max-1);
					int index_b = convert_triplet_to_index(_data, b, t, d_max-1);
					IloExpr expr(env);
					expr = Z[index_a] - Z[index_b];
					Precedence.add(IloRange(env, -1, expr, 0));
				}
		}
		
		// (14) 
		for(int b=0; b<nBlocks; b++)
			for(int t=0; t<t_max; t++)
				for(int d=1; d<d_max; d++){
					int index_l = convert_triplet_to_index(_data, b, t, d-1);
					int index_r = convert_triplet_to_index(_data, b, t, d);
					IloExpr expr(env);
					expr = Z[index_l] - Z[index_r];
					Precedence.add(IloRange(env, -1, expr, 0));

		}
		// (15)
		for(int b=0; b<nBlocks; b++)
			for(int t=1; t<t_max; t++){
				int index_l = convert_triplet_to_index(_data, b, t-1, d_max-1);
				int index_r = convert_triplet_to_index(_data, b, t, 0);
				IloExpr expr(env);
				expr = Z[index_l] - Z[index_r];
				Precedence.add(IloRange(env, -1, expr, 0));
		}
		
		IloRangeArray2 Resource(env);
		// (16) add constraints in the same order as lambda is constructed
		for(int t=0; t<t_max; t++){
			Resource.add(IloRangeArray(env));
			for( int r=0; r<r_max; r++){
				std::vector<double> coef(nZ, 0.0);
				std::vector<bool> ind(nZ, 0);/// could be used limits, instead
				for(int b=0; b<nBlocks; b++){
						int index;
						/// q_{1,r}z_{bt1}
						index = convert_triplet_to_index(_data, b, t, 0);
						if(ind[index]==0) ind[index] = 1; 
						coef[index] += (*blocks)[b].getRCoef(0,r);
						
						/// sum_{d=1}^{d_max} q_{dr}*(z_{btd}-z_{bt,d-1}) 
						for(int d=1; d<d_max; d++){
							/// z_{btd}
							index = convert_triplet_to_index(_data, b, t, d);
							if(ind[index]==0) ind[index] = 1; 
							coef[index] += (*blocks)[b].getRCoef(d,r);
							/// -z_{bt,d-1}
							index = convert_triplet_to_index(_data, b, t, d-1);
							if(ind[index]==0) ind[index] = 1; 
							coef[index] -= (*blocks)[b].getRCoef(d,r);
						}
						if(t>0){
							/// -q_{1r}*z_{b,t-1,d_max}
							index = convert_triplet_to_index(_data, b, t-1, d_max-1);
							if(ind[index]==0) ind[index] = 1; 
							coef[index] -= (*blocks)[b].getRCoef(0,r);
						}
				}
				
				IloExpr expr(env);
				for(int i=0; i<nZ; i++)
					if(ind[i]) expr += coef[i]*Z[i];
					
				Resource[t].add(IloRange(env, -IloInfinity, expr,_data.getLimit(r, t)));  ///only type L
			}
		}
		// Hz=h
		IloRangeArray Cuts(env);
		int n = (int) _Cpart.size();
		for(int i=0; i<n; i++){
			int m = (int) _Cpart[i].size();
			if(m>1) 
				for(int j=0; j<m-1; j++){
					IloExpr expr(env);
					expr = Z[_Cpart[i][j]] -Z[_Cpart[i][j+1]];
					Cuts.add(IloRange(env, 0, expr, 0));
				}
		}
		IloExpr exprObj(env);
		{
			std::vector<double> coef(nZ, 0.0);
			const double rate =_data.getDiscountRate();
			for(int b=0; b<nBlocks; b++){
				double cur_rate = 1.0;
				for(int t=0; t<t_max; t++){
					for(int d=0; d<d_max; d++){
						int index = convert_triplet_to_index(_data, b, t, d);
						coef[index] += (*blocks)[b].getProfit(d)*cur_rate; 
						if(d<d_max-1) coef[index] -= (*blocks)[b].getProfit(d+1)*cur_rate; 
						else
							if(t<t_max-1) coef[index] -= (*blocks)[b].getProfit(0)*cur_rate/(1+rate); 
					}
					cur_rate /= (1+rate);
				}
			}
			
			for(int i=0; i<nZ; i++) exprObj += coef[i]*Z[i];
		}
		IloObjective Obj=IloMaximize(env, exprObj);
		// cplex model 
		 IloModel restLP(env);
		 restLP.add(Precedence);
		 for(int t=0; t<t_max; t++)
			restLP.add(Resource[t]);
		 restLP.add(Cuts);
		 restLP.add(Obj);
		 
		 IloCplex rlp_cpx(restLP);
		 rlp_cpx.setOut(env.getNullStream());
		 rlp_cpx.setParam(IloCplex::ClockType,1);
		 rlp_cpx.solve();
		 if(rlp_cpx.getStatus() == IloAlgorithm::Optimal) _status = 1;
		 else _status = 0; 
		 if(_status){
			IloNumArray sol_rlp(env);
			rlp_cpx.getValues(sol_rlp, Z);
			for(int i=0; i<nZ; i++)	_sol[i] = sol_rlp[i];
			_objVal = (double) rlp_cpx.getObjValue();
			for(int t=0; t<t_max; t++){
				rlp_cpx.getDuals(sol_rlp, Resource[t]);
				for( int r=0; r<r_max; r++)	_pi[t][r]= sol_rlp[r];
			}
		}
	// free memory
		for(int t=0; t<t_max; t++){
			Resource[t].endElements();
			Resource[t].end();
		}
		Resource.end(); 
	}
	catch (IloException & ex) { cerr << "Error: " << ex << endl; }
	catch (...) { cerr << "Unknown Error Occured" << endl;}
	timer.stop();
	_t_cpu = timer.getTime();
	env.end();
}
bool ResLP::getStatus(){return _status;}

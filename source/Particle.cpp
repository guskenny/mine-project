#include "Particle.h"
#include "VolVolume.hpp"
//#include "Random.h"
#include "MaxClosureFactory.h"
#include "SinglePSolver.h"
#include "UpitSolver.h"
#include <queue>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp> // define boost::tie

using namespace LaPSO;

ResourceUse::ResourceUse(const CumulativeModel &prob)
{
    resize(prob.graph.getNumNodes());
    for(size_t r=0;r<prob.resLim.size();++r){
	CumulativeModel::Vertex v;
	double q;
	BOOST_FOREACH(boost::tie(v,q),prob.res[r])
	    (*this)[v].push_back( std::make_pair(r,q) );
    }
}

bool ResourceUse::canFit(CumulativeModel::Vertex v,const std::vector<double> &residual) const
{
    int r; double q;
    BOOST_FOREACH(boost::tie(r,q),(*this)[v]){
		if(residual[r] < q ) return false;
    }
    return true;
}
// update residual when vertex v is used
void ResourceUse::update(CumulativeModel::Vertex v,std::vector<double> &residual) const
{
    int r; double q;    
    BOOST_FOREACH(boost::tie(r,q),(*this)[v]){
	residual[r] -= q;
    }
}
// undo update residual when vertex v is used
void ResourceUse::undoUpdate(CumulativeModel::Vertex v,std::vector<double> &residual) const
{
    int r; double q;    
    BOOST_FOREACH(boost::tie(r,q),(*this)[v]){
	residual[r] += q;
    }
}


// LaPSO solves max_dual min_x -profit * x + dual (resLim - resUse(x)) + perturbation
LaPSO::Status PCPSPpsoHooks::reducedCost(const LaPSO::Particle &p_, LaPSO::DblVec &redCost)
{
    const PCPSPparticle &p(static_cast<const PCPSPparticle &>(p_));
    redCost = prob.getProfit();
    redCost *= -1.0/SCALE;
    for(size_t r=0;r<prob.resLim.size();++r){
	CumulativeModel::Vertex v; double q;
	BOOST_FOREACH(boost::tie(v,q),prob.res[r]){
	    redCost[v] -= q*p.dual[r]/SCALE;
	}
    }
    return LaPSO::OK;
}


// LaPSO solves max_dual min_x -profit * x + dual (resLim - resUse(x))
// where dual <= 0 
LaPSO::Status PCPSPpsoHooks::solveSubproblem(Particle &p_)
{
    PCPSPparticle &p(static_cast<PCPSPparticle &>(p_));
    p.rc *= -SCALE;		// MaxClosure maximises rc is for minimisation
    p.mc->setProfit(p.rc);
    p.mc->solve(); ++nMCsolves;
    //std::vector<int> sol(prob.graph.getNumNodes());
    p.mc->getClosure(p.x);
    // p.lb = relaxed bound we will invert sign at the end to get -ve profit
    p.lb = p.mc->calcProfit(p.rc,p.x); // (profit + dual*resUse)*x
    // add lagrangian constant:
    for(size_t r=0;r<prob.resLim.size();++r)
	p.lb -= prob.resLim[r]*p.dual[r];
    // upper bound: inverted sign (we are minimising negative objective
    p.ub = p.mc->calcProfit(prob.getProfit(),p.x); // original profit
    p.lb /= -SCALE;
    p.ub /= -SCALE;
    p.rc *= -1.0/SCALE;		// return to original value
   
    // subgradient / violation calculation
    p.isFeasible = true;
    for(size_t r=0;r<prob.resLim.size();++r){
	p.viol[r] = prob.resLim[r];
	CumulativeModel::Vertex v;
	double q;
	BOOST_FOREACH(boost::tie(v,q),prob.res[r])
	    if(p.x[v]==1) p.viol[r] -= q;
	if(p.viol[r] > 1e-7*fabs(prob.resLim[r]))
	    p.isFeasible = false;
	p.viol[r] /= SCALE;
    }
    std::cout<<"\tLag/10^6=" << -p.lb << ",  infeasCost *1e-6=" << -p.ub
	     << ", lag=" << p.dual.min() << " to " << p.dual.max()
	     << "\n\tViol = " << p.viol.min() << " to " << p.viol.max() 
	     << std::endl; 
    
    // adjust lower bound for perturbation:
    // Since reduced cost is perturbed we may have underestimated profit
    // we know
    //  (profit-perturb)*optimal <= (profit-perturb) * x
    // hence   profit * optimal <= (profit-perturb)*x + perturb*optimal
    //                          <= (profit-perturb)*x + max_closure(perturb)
    // However too expensive to calculate max_closure(perturb). Can be
    // bounded by sum( max(0,perturb) )
    // rc = reduced cost = (-profit+perturb)/SCALE
    // lb = -profit*optimal/SCALE >= rc * x - sum( max(0,perturb) )
    for(size_t i=0;i< p.perturb.size();++i){
	p.lb -= std::max(0.0,p.perturb[i]);
    }
    return LaPSO::OK;	// could return ABORT to stop early
} 

LaPSO::Status PCPSPpsoHooks::fixConstraint(const int res,const Particle &p_,SparseVec &feas)
{
    return NONE; // should really return a maximal set of vertices that form
		 // closure and respect single resource constraint but that's
		 // too hard for now
}
LaPSO::Status PCPSPpsoHooks::heuristics(Particle &p_)
{
    PCPSPparticle &p(static_cast<PCPSPparticle &>(p_));
    // repeated greedy construction based on current state
    DblVec resAvail(prob.resLim);
    std::vector<int> blockCnt(prob.graph.getNumNodes());
    // avail is available blocks sorted largest to smallest
    std::priority_queue<std::pair<double,CumulativeModel::Vertex> > avail;
    IntVec soln(prob.graph.getNumNodes(),0);
    DblVec priority(prob.graph.getNumNodes());
    IntVec x = p.x; // current subproblem solution
    p.ub = 1e99;
    p.isFeasible=false;
    int bestTrial=-1,trial;
	qol::CpuTimer timer;
    for(trial=0;trial < 16;++ trial ) { // was 8
    	// try constructing a better solution guided by a priorty

		//---------------- set priority ---------------------------------
		double perturbSize= std::min( // want a smallish number
				0.8,std::max(1e-6,p.rand()*p.rand()));
		switch(trial %4 ){
		case 0: priority=p.perturb;	// 'random' perturbation
			priority.negate(); // perturbation +ve ==> penalise inclusion
			perturbSize *= 1e-2; // already perturbed, only break ties
			break;
		case 1: priority=p.rc;
			priority.negate(); // reduced cost is for minimisation
			perturbSize *= 0.1;
			break;
		case 2:
			priority = x; break; // integer soln
		case 3:
			priority = prob.getProfit();
			break;
		}
		// now add a random perturbation
		double minPri=priority.min(),maxPri=priority.max();
		for(size_t i=0;i<priority.size();++i){
			priority[i] += p.rand()*perturbSize*(maxPri - minPri);
		}
		const int excludeZero= (trial/3)%2; // don't include v if x == zero

		//----------------- reset -------------------------
		resAvail = prob.resLim;
		soln=0;
		// avail.clear(); // always empty at this point anyway
		for(int v=0;v<prob.graph.getNumNodes();++v){
			blockCnt[v] = prob.graph.getNode(v)->getInDegree();
			if(blockCnt[v] == 0 &&  // accessible
			   x[v] >= excludeZero) // only consider blocks in current solution
				avail.push( std::make_pair(priority[v],v) );
		}
		long tried=0,accepted=0;
		double obj = 0;
		//----------- add accessible blocks if feasible -----------
		while( ! avail.empty() ){
			CumulativeModel::Vertex vv = avail.top().second;
			avail.pop();		// remove top element
			if(soln[vv] == 1) continue;	// already included
			for(CumulativeModel::Vertex v=vv;
				v!= InvalidVertex;
				v=prob.getSuccVertex(v)){
				if(blockCnt[v] != 0) break;
				++tried;
				if( resUse.canFit(v,resAvail) ){
					soln[v] = 1;
					vv = v;
					++accepted;
					obj += prob.getProfit(v);
					resUse.update(v,resAvail); // remove resource usage
					const Node *nd = prob.graph.getNode(v);
					for(int i=0;i<nd->getOutDegree(); ++i){
						int w = nd->getOutArc(i)->getTgtID();
						--blockCnt[w]; 
						if(blockCnt[w] == 0 && x[w] >= excludeZero
						   && w != prob.getSuccVertex(v)) // will try anyway
							avail.push(std::make_pair(priority[w],w));
					}
				}else
					break;
			} // end loop over vertices for same block
			while(soln[vv] == 1 && prob.graph.getNode(vv)->getOutDegree() <= 1
				  && prob.getProfit(vv) < 0 ) { // don't include this vertex
				soln[vv] = 0;
				obj -= prob.getProfit(vv);
				--accepted;
				if( prob.graph.getNode(vv)->getOutDegree() == 1)
					++blockCnt[prob.graph.getNode(vv)->getOutArc(0)->getTgtID()];
				resUse.undoUpdate(vv,resAvail);
				vv = prob.getPredVertex(vv);
				if(vv == InvalidVertex) break;
			}
		}
		obj /= SCALE;
		if( -obj < p.ub){
			bestTrial = trial;
			p.ub = -obj;// upper bound of -ve objective function
			p.x = soln;
			p.isFeasible = true;    
		}
		//if(trial%4==0)
		//	printf("\t%d: heuristic %g with %ld/%ld accepted, priority %g to %g\n",
		//		   trial,obj,tried,accepted,minPri,maxPri);
    } // end loop over trials
    printf("Heuristic best = %g (global %g) found at trial %d/%d in %.2f sec\n",
		   -p.ub,bestCost/SCALE,bestTrial,trial,timer.elapsedSeconds());
    return OK; 			// found a solution - possibly a rubbish one
    
}
LaPSO::Status PCPSPpsoHooks::updateBest(Particle &p)
{
    bestCost = -p.ub*SCALE;
	bestSoln = p.x;
	printf("#### New best %.8g found after %.2f sec wall %.2f sec CPU\n",
		   bestCost,timer.elapsedWallTime(),timer.elapsedSeconds());
    return OK;
}

// main function for Lagrangian particle swarm
LaPSO::Problem *runParticleMain(CumulativeModel &prob,int argc,const char **argv)
{
    const int nVert = prob.graph.getNumNodes(),nRes=prob.res.size();
	LaPSO::Problem *lapso=new LaPSO::Problem(nVert,nRes);
    LaPSO::Problem &solver=*lapso; // num variables/constraints
    PCPSPpsoHooks pcpsp(prob);
    //-------- set default parameter values
    solver.param.absGap = 0.999;
    solver.param.printFreq = 1;	
    solver.param.printLevel = 1;
    solver.param.maxIter = 200000;
    solver.param.subgradFactor = 0.2; // seems OK
    solver.param.heurFreq=1;	// heuristic is fast compared to lagrangian
    bool useVol = false;
    bool useUpit = false; // flag to use UpitSolver
    MaxClosureFactory mcfactory;    
    for(int i=0;i<argc;++i){
	if(argv[i][0] == '-' && argv[i][1] != '-')
	    switch(argv[i][1]){
		case 't':
		    solver.param.maxWallTime=atof(argv[i+1]); break;
		case 'V': useVol = true; break;
		case 'U': useUpit = true; break;
		case 'n': case 'b': case 'e': case 'p':
		case 'N': case 'B': case 'E': case 'P':
		    mcfactory.setOption(argv[i][1]);
		    break;
	    }
    }
    // override defaults with command line arguments
    solver.param.parse(argc,argv);
    double maxCost = *std::max_element(prob.getProfit().begin(),prob.getProfit().end());
    maxCost /= SCALE;
    pcpsp.bestCost = solver.best.ub = 0;
    solver.best.lb = -1e20; //nothing interesting yet
    solver.dualUB = 0;	    // lagrange multipliers must be <= 0
    Uniform rand;
    if(solver.param.randomSeed == 0) rand.seedTime();
    else rand.seed(solver.param.randomSeed);

    DblVec initLag(2,0); double norm=0;
    for(size_t r=0;r<initLag.size();++r){
		CumulativeModel::Vertex v;
		double q;
		initLag[r] = +prob.resLim[r];
		BOOST_FOREACH(boost::tie(v,q),prob.res[r])
			initLag[r] -= q;
		norm+= initLag[r]*initLag[r];
    }

    std::map<int,int> fixed;

    // run UpitSolver to find any fixable variables
    if(useUpit){
      UpitSolver uSolve = UpitSolver(prob);
      //SinglePSolver sp_solver = SinglePSolver(prob, argc, argv);
      uSolve.solve();
      uSolve.getFixed(fixed);
      std::cout << "Found " << fixed.size() << " variables to fix" << std::endl;
    }

    norm = sqrt(norm);
    for(size_t r=0;r<initLag.size();++r) initLag[r]= -fabs(initLag[r])/norm;
    if(useVol){
		PCPSPparticle *p = new PCPSPparticle(nVert,nRes,mcfactory(prob.graph,fixed));
		p->rand.seed(std::numeric_limits<uint32_t>::max()*rand());
	    for(int j=0;j<solver.dsize;++j)
			p->dual[j] = (j < initLag.size()) ? rand(0.5,1.0)*initLag[j] :
				rand(0.7,1.0)*p->dual[j-initLag.size()];
		LaPSO::VOL_LaPSO_adaptor vol(solver,pcpsp,p);
		vol.prob.parm.lambdainit = 10; //solver.param.globalFactor;
		vol.prob.parm.greentestinvl = 3;
		vol.prob.parm.redtestinvl = 5;
		std::cout << "Using volume algorithm with random initial point\n";
		vol.solve(true);
		solver.best = vol.best;
    }else{
		solver.swarm.resize(solver.param.nParticles,0);
#       pragma omp parallel for 
		for(int i=0;i<solver.param.nParticles;++i)
			solver.swarm[i] = new PCPSPparticle(nVert,nRes,mcfactory(prob.graph,fixed));
		std::cout << "Created " << solver.param.nParticles << " particles"
				  << " in " << pcpsp.timer.elapsedWallTime() << " sec\n";
		for(int i=0;i<solver.param.nParticles;++i){ // serial rand seed init
			PCPSPparticle *p=(PCPSPparticle *)solver.swarm[i];
			p->rand.seed(std::numeric_limits<uint32_t>::max()*rand());
		}
#       pragma omp parallel for 
		for(int i=0;i<solver.param.nParticles;++i){ // init dual values
			PCPSPparticle *p=(PCPSPparticle *)solver.swarm[i];
			for(int j=0;j<solver.dsize;++j)
				p->dual[j] = (j < initLag.size()) ?
					p->rand(0.5,1.0)*initLag[j] :
					p->rand(0.7,1.0)*p->dual[j-initLag.size()];
		}
	    //solver.swarm.push_back(p);
		std::cout << "set up solver with " << solver.swarm.size() << " particles\n";
		solver.solve(pcpsp);
	}
    std::cout << "Completed in " << pcpsp.timer.elapsedWallTime() << " sec, "
			  << pcpsp.timer.elapsedSeconds() << " CPU sec "
			  << pcpsp.nMCsolves << " MaxClosure calls "
			  << 100.0*pcpsp.timer.elapsedSeconds()/(
				  solver.param.nCPU*pcpsp.timer.elapsedWallTime())
			  << "% utilisation\n";
    if(solver.best.isFeasible){
		std::cout << "best solution = " << -SCALE*solver.best.ub
				  << " <= " << -SCALE*solver.best.lb << std::endl;
		/* should dump out solution here somewhere
		   std::cout << "F:";
		   for(int i=0;i<dcmst.n-1;++i) printf("%3d",dcmst.mst[i].i);
		   std::cout <<"\nT:";
		   for(int i=0;i<dcmst.n-1;++i) printf("%3d",dcmst.mst[i].j);
		   std::cout << std::endl;
		*/
    }else
		std::cout << "no feasible solution found\n";
	return lapso;	
} // end runParticleMain()

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/

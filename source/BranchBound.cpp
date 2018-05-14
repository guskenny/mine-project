#include "BranchBound.h"
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp> // define boost::tie
#include "MaxClosure_Base.h"
#include "MaxClosure_BoostMaxFlow_BK.h"
#include "MaxClosure_BoostMaxFlow_EK.h"
#include "MaxClosure_BoostMaxFlow_PR.h"
#include "MaxClosure_NetworkFlow.h"
#include "graph.h"
#include "QOL/CpuTimer.h"
#include "MaxClosureFactory.h"
#include "Particle.h"		// LaPSO & Volume algorithms
#include "bz.h"


void NodeSolver::LagrangianTest()
{
    std::map<int,int> fixed;
    prob.BN_info_to_map(node,fixed);    
    MaxClosureFactory mcfactory(maxClosureMethod);
    bool debug=true;
#   ifdef NDEBUG
    debug = false;
#   endif    
    qol::CpuTimer tmp;
    MaxClosure_Base *mc=mcfactory(prob.graph,fixed);
    std::cout << "MaxClosure " << mcfactory.getOption() << " setup in "
	      << tmp.elapsedSeconds() << " sec\n"; std::cout.flush();
    std::vector<double> lag(prob.getnConstraints(),0.0);
    std::vector<double> rc(prob.graph.getNumNodes());
    std::vector<double> viol(prob.getnConstraints());
    std::vector<int> soln(rc.size(),1); // initial solution - everything mined immediately
    qol::CpuTimer timer; int iter,maxiter;
    if(timeLimit <= 0) {
	maxiter = 10;
	timeLimit = 1e20;
    }else
	maxiter= 99999;
    timer.setTimeLimit(timeLimit);
    double UB=1e99;
    for(iter=0;iter<maxiter && ! timer.timeLimitReached();++iter){
	std::cout << "Lagrangian test iteration " << iter << " --------------\n";
	{ // really crude subgradient iteration
	    double step = 1.0/(1+0.5*iter),maxSlack=0,maxViol=0,norm=0;
	    for(size_t r=0;r<prob.resLim.size();++r){
		viol[r] = 0;
		CumulativeModel::Vertex v;
		double q;
		BOOST_FOREACH(boost::tie(v,q),prob.res[r])
		    if(soln[v]==1) viol[r] += q;
		if(debug && viol[r] > 0)
		    std::cout << (boost::format("\tC%02d(%03d): %.2f <= %.2f\n")%r
				  %(int)prob.res[r].size()%viol[r]%prob.resLim[r]);	
		viol[r] -= prob.resLim[r];
		if(viol[r] > maxViol) maxViol = viol[r];
		if(viol[r] <-maxSlack) maxSlack= -viol[r];
		norm += viol[r]*viol[r];
	    }
	    step /= sqrt(norm);
	    std::cout <<"\tMax violation = " << maxViol << std::endl;
	    std::cout <<"\tMax slack     = " << maxSlack << std::endl;
	    std::cout <<"\tstep size     = " << step << std::endl;
	    for(size_t r=0;r<prob.resLim.size();++r)
		lag[r] = std::max(0.0,lag[r]+step*viol[r]);
	    // In the first iteration use the violation of the first constraint
	    // (only one violated) to guide initial values for all others
	    if(iter==0) for(size_t r=1;r<lag.size();++r) lag[r]=lag[r-1]/(1+prob.getDiscountRate());
	}
	double lagConst = prob.calcRedCost(lag,rc);
	{ double maxRC=-1e99,minRC=1e99;
	    for(size_t i=0;i<rc.size();++i){
		maxRC=std::max(maxRC,rc[i]);
		minRC=std::min(minRC,rc[i]);
	    }
	    std::cout << boost::format("\tRed cost range=%g .. %g\n")%minRC%maxRC;
	}
	std::cout.flush();	  
	double start=timer.elapsedSeconds();//,wall=timer.elapsedWallTime();
	if( iter > 0 && start/iter + start > timeLimit){
	    std::cout << "Insufficient time for another iteration\n";
	    break;
	}
	mc->setProfit(rc);
	std::cout << "\tSolve status: " << mc->solve() << std::endl;
	std::cout << "\tGet closure:  " << mc->getClosure(soln) << std::endl;
	double obj=mc->calcProfit(rc,soln);
	int cnt=0;
	for(size_t v=0;v<soln.size();++v) if(soln[v]){ ++cnt;}
	std::cout << (boost::format("\tLag. UB = %.4e + %.4e = %g")
		      %lagConst%obj%(obj+lagConst))
		  << (boost::format(" in %.2f sec CPU\n")%
		      (timer.elapsedSeconds()-start)) //%(timer.elapsedWallTime()-wall))
		  << (boost::format("\t%.2f total sec. %d vertices included\n")
		      %timer.elapsedSeconds()%cnt);
	UB = std::min(UB,obj+lagConst);
	if(debug) // isClosure() takes several seconds on large instances
	    std::cout << "\tClosure validity: " << mc->isClosure(soln)<< std::endl;
	// double check profit
	std::vector<double> fsoln(soln.size());
	for(size_t i=0;i<soln.size();++i) fsoln[i] = soln[i];
	double ob=0.0;
	std::cout << "\tPeriod: no. blocks processed/mined\n\t";
	for(int p=0;p<prob.getNPeriod();++p){
	    int blockCnt=0,procCnt=0; // processed = destination 1
	    for(int b=0;b<prob.getNBlock();++b){
		double mined=prob.blockMined(fsoln,b,p);
		double total=0;
		for(int d=0;d<prob.getnDestination();++d){
		    const double dest =prob.blockProcessed(fsoln,b,d,p);
		    if(dest != 0){
			ob += prob.getProfit(b,d,p) * dest;
			total += dest;
			if(d==1) ++procCnt;
		    }
		}
		if(fabs(total - mined)> 1e-5){
		    std::cerr <<"ERROR: destination total doesn't match "
			      << boost::format("b=%d,t=%d - %f vs %f")%b%p%mined%total
			      << std::endl;
		}
		blockCnt += int(mined+0.1);
	    }
	    if(blockCnt > 0)
		std::cout << boost::format("%d:%d/%d ")%p%procCnt%blockCnt;
	}
	std::cout << "\n\tPrimal (non-feas.) profit is " << ob << std::endl;
	std::cout.flush();
    }
    std::cout << "Completed " << iter << " iterations in "
	      << (boost::format("%.2f / %.2f sec CPU/wall, best upper bound = %g\n")
		  %timer.elapsedSeconds()%timer.elapsedWallTime()%UB);
} // end LagrangianTest
  


void NodeSolver::solve() // updats sol_int, and possibly sol_real
{
    MaxClosureFactory mcfactory(maxClosureMethod);// NEED TO ADD FIXED HERE
    switch(runMethod){
	case 'T': LagrangianTest(); break;
	case 'V': case 'L':
	    runParticleMain(prob,arg.size(),&arg[0]);
	    break;
	case 'Z': {
	    BZ bz(prob,mcfactory,timeLimit,&node); 
	    bz.solve();
	    break;
	}
	default:
	    std::cerr << "WARNING: unkown run method: " << runMethod << std::endl;
	    runParticleMain(prob,arg.size(),&arg[0]);
    }

}

void NodeSolver::parseArgs(int argc,const char **argv)
{
    int opt;
    runMethod='T';		// test
    int narg=argc-1;
    arg.resize(narg);
    for(int a=0;a<narg;++a) arg[a] = argv[a];
    // skipp any long arguments with --
    while(narg > 0 && argv[narg-1][0]=='-' && argv[narg-1][1]=='-') --narg;
    //std::cout << "Parsing " << narg << " arguments\n";
    while ((opt = getopt(narg, (char *const*)argv, "hnebpw:NEBPt:VTLZ ")) != -1) {
	switch ((char)opt) {
	case 'n': case 'b': case 'e': case 'p':
	    case 'N': case 'B': case 'E': case 'P':
	    maxClosureMethod = (char)opt;
	    break;
	case 't': timeLimit = atof(optarg); break;
	case 'V': case 'L': case 'Z':
	    runMethod = (char)opt; break;
	default:
	    continue;		// ignore additional options (for runParticleMain)
	}
    }

}


int BranchBoundTree::callSolver(MineProblem &tree)
{
    const double relGap = 0.01;
    if( ! solver ) solver = new SolverInterface(*data);
    best = *tree.current; // store best solution
    int status=0,iter=0;
    while(tree.openSubProblems() > tree.solvedSubProblems() && tree.current != 0){
	std::cout << "Solving branch & bound node " << iter++
		  << " open/total = " << tree.openSubProblems() << "/" << tree.solvedSubProblems()
		  << std::endl;
	qol::CpuTimer timer;
	if(alg_type==99){
	    NodeSolver node(*data,best);
	    node.LagrangianTest();
	}else
	  status = solver->callSolver(tree,0,1,alg_type,-1,-1); // call BZ algorithm
	if(tree.current->sol_int.obj > best.sol_int.obj)
	    best.sol_int = tree.current->sol_int;	    
	double LB=best.sol_int.obj, UB=tree.current->sol_real.obj;
	std::cout << "Tree solver LB="<<LB  << ",  UB=" << UB
		  << ", gap=" << 100*(UB-LB)/LB << "%  "
		  << tree.current->sol_int.numMined() << " blocks mined\n";
	std::cout << boost::format("Solution time = %.2f CPU, %.2f sec elapsed time\n"
	    )%timer.elapsedSeconds()%timer.elapsedWallTime();
	tree.setSolved();
	if(tree.current->sol_real.obj > (1.0+relGap)*best.sol_int.obj){
	    // split problem here
	    ;
	}
	tree.next();
    }
    return status;
}



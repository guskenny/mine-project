#ifndef __BRANCH_BOUND_H__
#define __BRANCH_BOUND_H__

#include "CumulativeModel.h"
#include "../parallel/MineProblem.h" // defines ProblemData
#include "solver_interface.h"
#include <list>

class NodeSolver : public ProblemData {
    // NodeSolver includes node, sol_int, sol_real
public:
    CumulativeModel &prob;
    /// Possible runMethods include
    /// A : ACO-MIP (not yet done)
    /// L : LaPSO
    /// T : Simple test lagrangian method
    /// V : Volume method
    /// Z : Bienstock-zuckerberg
    char runMethod, maxClosureMethod;
    double timeLimit; // wall-clock time limit
    int maxThreads;   // defaults to 1
    std::vector<const char *> arg;

    NodeSolver(CumulativeModel &_prob,ProblemData &_branch,
	       int argc=0,const char **argv=0) // optional arguments
	: prob(_prob), ProblemData(_branch),timeLimit(-1),maxThreads(1),arg(argc)
	{ for(int i=0;i<argc;++i) arg[i] = argv[i];}
    virtual void solve(); // updats sol_int, and possibly sol_real

    void LagrangianTest();
    void parseArgs(int argc,const char **argv);
};


class BranchBoundTree : public SolverInterface {
public:
    SolverInterface *solver;			  // subproblem solver
    ProblemData best;
    int alg_type;
    BranchBoundTree(CumulativeModel &data) : SolverInterface(data),solver(0) {}
    virtual ~BranchBoundTree() { if(solver) delete solver; }
    virtual int callSolver(MineProblem &problem); // executes branch & bound
    
}; // end BranchBoundTree
    

// main function for a branch and bound algorithm
// can pass commandline arguments in to control what method to use for
// improving the branch and bound performance
/* inline void runBranchBound(CumulativeModel &prob,
		    ProblemData &branch, // input: current state, output: best soln found
		    int argc=0,const char **argv=0)
{ BranchBoundTree tree(prob,branch,argc,argv);
  tree.solve();
  } */

#endif

//---------------------------------------------------------------------------------
// Main program
// Note: Most of the material from here has moved to BranchBound.cpp
// (NodeSolver::solve() method implements solution for a single branch & bound node)

#include <stdlib.h>
#include "CumulativeModel.h"
#include "../parallel/MineProblem.h" // defines ProblemData 
#include "BranchBound.h"

void usage(const char **argv)
{
    std::cerr << argv[0] << " [options] Path\n"
	      << "Solve PCPSP (mine planning optimisation) using a cummulative variable formulation\n"
              << "Path    directory containing data files (uncompressed)\n"
              << "-C      CPLEX solver (not implemented yet)\n"
              << "-G      Gurobi solver\n"
              << "-X      eXact branch and bound tree method\n"
              << "-x n    limit number of branch and bound nodes\n"
	      << "-M      MIP-ACO"
              << "-V      Volume method \n"
              << "-L      Lagrangian particle swarm optimisation\n"
              << "-T      Lagrangian test (very simplistic optimisation) \n"
              << "-n/N    network flow max closure solver / with pre-processing\n"
              << "-b/B    Boykov-Kolmogorov max flow algorithm for max closure\n"
              << "-e/E    Edmonds-Karp max flow algorithm\n"
              << "-p/P    Push-relable max flow algorithm\n"
              << "-w fn   Write cumulative model to file 'fn'\n"
              << "-Z      run bienstock-zuckerberg LP method\n"
	      << std::endl;
}


int main(int argc,const char *argv[]){

    if(argc <= 1){
	usage(argv);
	return 1;			// failed
    }
    CumulativeModel prob(argv[argc-1]);
    std::cout << "Loaded " << prob.getName()
	      << (boost::format(" with %d blocks, %d destinations, %d periods,")%
		  prob.getNBlock()%prob.getnDestination()%prob.getNPeriod())
	      << (boost::format(" %d resources & %d arcs")
		  %prob.getnResources()%prob.getnArcs())
	      <<std::endl;
    std::cout << (boost::format("Cumulative problem has %d nodes, %d arcs, %d constraints")
		  %prob.graph.getNumNodes()%prob.graph.getNumArcs()
		  %prob.getnConstraints()) << std::endl;
    int opt, narg=argc-1;
    while(narg > 0 && argv[narg-1][0]=='-' && argv[narg-1][1]=='-') --narg;
    int alg_type=2;
    while ((opt = getopt(narg, (char *const*)argv, "hnebpw:GMNEBPt:VTLZx:X ")) != -1) {
	switch ((char)opt) {
	    /* -- these options only parsed in Solve
	      case 'n': case 'b': case 'e': case 'p':
	case 'N': case 'B': case 'E': case 'P':
	    maxClosureMethod = (char)opt;
	    break;
	    case 't': timelimit = atof(optarg); break;
	case 'V': case 'L': case 'Z':
	runMethod = (char)opt; break; */
	case 'M': alg_type=5; break;
	case 'V': alg_type=6; break;
	case 'L': alg_type=7; break;

	case 'G': alg_type=8; break;
	case 'Z': alg_type=2; break;
	case 'T': alg_type=99; break; // test
	case 'w':			// write to file
	    prob.dump(optarg);
	    std::cout << "Wrote " << optarg << std::endl;
	    break;
	case 'h':		// help
	    usage(argv);
	    return 2;
	default:
	    break;		// ignore additional options (for runParticleMain)
	}
    }
    //BranchNode_info branch(prob.getNBlock(),prob.getNPeriod(),prob.getnDestination(),prob.getnResources());
    //ProblemData branch(prob.getNBlock(),prob.getNPeriod(),prob.getnDestination(),prob.getnResources());
    
    //for(int b=0;b<500;++b) // do first n blocks in period 0
    //  branch.time[b][1]=1;

    //BranchBoundTree tree(prob,branch,argc,argv);
    BranchBoundTree tree(prob);
    tree.alg_type = alg_type;
    MineProblem worklist(argv[argc-1],0,1); // create empty tree
    worklist.next();
    worklist.current->init(prob.getNBlock(),prob.getNPeriod(),prob.getnDestination(),prob.getnResources());
	/*
    ProblemData *root = new ProblemData(
	prob.getNBlock(),prob.getNPeriod(),prob.getnDestination(),prob.getnResources());
    std::cout <<"x size=" << root->sol_int.x.size() << std::endl;
    root->assignedProc=0;
    worklist.addWork(root );
    worklist.current = root;
	*/
    std::cout <<"x size=" << worklist.current->sol_int.x.size() << std::endl;
    tree.callSolver(worklist);
    std::cout << "done\n";
    return 0;  
}

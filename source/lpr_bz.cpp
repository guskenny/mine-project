/* LP Relaxation: Open Pit Mining(PCPSP)
 * author: Davaa Baatar
 * date: 19.06.2016
*/
#include <iomanip>
#include <stdlib.h>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp> // define boost::tie
#include "CumulativeModel.h"
#include "MaxClosureFactory.h"
#include "../QOL/CpuTimer.h"
#include "util.h"
#include "rlp.h"
#include "bz.h"
#include "../LagrangianHeuristic/solver.h"
#include "../LagrangianHeuristic/solver_functions.h"
#include "../LagrangianHeuristic/pheromones.h"
#include "../LagrangianHeuristic/ACO_solution.h"

int main(int argc, char** argv){
	if(argc <= 1){
		usage(argv);
		return 1;			// failed
	}
	
	// read data and construct the graph
	CumulativeModel prob(argv[argc-1]);
	std::cout << "Loaded " << prob.getName()
			  << (boost::format(" with %d blocks, %d destinations, %d periods,")
					%prob.getNBlock()%prob.getnDestination()%prob.getNPeriod())
			  << (boost::format(" %d resources & %d arcs")
					%prob.getnResources()%prob.getnArcs())
			  << std::endl;
	std::cout << (boost::format("Cumulative problem has %d nodes, %d arcs, %d constraints")
					%prob.graph.getNumNodes()%prob.graph.getNumArcs()
					%prob.getnConstraints()) 
			  << std::endl;
	
	/* Initialize */
	
	qol::CpuTimer timer;
	double cpu_time = timer.elapsedSeconds();
	MaxClosureFactory mcfactory;
	int opt;
	double timelimit = -1;
	
	while ((opt = getopt(argc-1, (char *const*)argv, "nebpw:NEBPt:")) != -1) {
		std::cout<<"Max Closure Solver: ";
		switch ((char)opt) {
		  case 'n': case 'N': 
				std::cout<<"Network Simplex Algorithm "<< std::endl; 
				mcfactory.setOption((char)opt);
				break;
		  case 'b': case 'B': 
				std::cout<<"Boykov-Kolmogorov max flow algorithm "<< std::endl; 
				mcfactory.setOption((char)opt);
				break;
		  case 'e': case 'E': 
				std::cout<<"Edmonds-Karp max flow algorithm "<< std::endl; 
				mcfactory.setOption((char)opt);
				break;
		  case 'p': case 'P': 
				std::cout<<"Push-relable max flow algorithm "<< std::endl; 
				mcfactory.setOption((char)opt);
				break;
		  case 't': timelimit = atof(optarg); break;
		  case 'w':			// write to file
			  prob.dump(optarg);
			  std::cout << "Wrote " << optarg << std::endl;
			  break;
		  default:
			  std::cerr << "Unknown option '" << (char) opt << "'\n";
			  usage(argv);
			  return 2;
		}
	}
	
	BZ lpr_bz(prob, mcfactory, timelimit);
	int status = 0;
	
	status = lpr_bz.solve();
	
	cpu_time = timer.elapsedSeconds() - cpu_time;
	
	std::cout<<(boost::format("\n%-36s %s") %"Status" %( (status) ? "Optimal " : "Feasible/NotSolved" ))<<std::endl;
	std::cout<<(boost::format("%-36s %-8.2f sec.") %"Total CPU time" %cpu_time)<<std::endl;
	if(status){
			const int nBlocks = prob.getNBlock();
			const int t_max = prob.getNPeriod();
			const int d_max = prob.getnDestination(); 
		
      Sol_Real sol_xy(nBlocks, t_max, d_max);  
      std::vector<std::vector<double> > &X = sol_xy.x;
			std::vector< std::vector< std::vector<double> > > &Y = sol_xy.y; 
			std::vector<double> & sol_z = lpr_bz.getSolution();
			int nPeriods = 0;
			//block -> time -> destination -> value
			convert_Z_to_XY(prob, sol_z, X, Y, nPeriods);/// X{bt} usage Y{b,t,d} destribution

			// run some heursitics
			// Create a MaxClosure object
			cout << "Running heuristics ... " << endl;

			Daten *data = new Daten(argv[argc-1],'p');
			SolverFunctions *funcs = new SolverFunctions(*data,1000); // diff_amt - fairly arbitrary
			Solver *alg_obj = new Solver(*data, *funcs);
			Pheromones *p = new Pheromones(data, funcs->profit);
			ACO_Solution * lr_sol = new ACO_Solution(data,funcs->profit,funcs->block_dest,funcs->succ);
			std::vector<bool> negBlocks(nBlocks);  // should do something more sophisticated here
			for(int b=0;b<nBlocks;++b) negBlocks[b] = (data->getBlock(b).getProfit(0) < 0);
			double obj = lr_sol->setFracBlocks(X,negBlocks);
			cout << "Best solution: " << lr_sol->getObj() << endl;

			ACO_Solution* a = alg_obj->ACO_run(p, lr_sol, 1000, 20, 0.2, 0.1, true, 1000);

			cout << "best ACO solution: " << a->getObj() << endl;
			
      sol_xy.obj = evaluateObjValue(*data, sol_xy);

			//delete a;
			delete lr_sol;
			delete p;
			delete funcs;
			delete alg_obj;
			delete data;
			// display info, save solution
			double obj_val = sol_xy.obj;
			double ub = lpr_bz.getUB();
			int nIter = lpr_bz.get_nIter();
			double total_MC_time = lpr_bz.getCpuTime_MC();
			double total_RLP_time = lpr_bz.getCpuTime_RLP();
			std::cout<<(boost::format("   %-33s %8.2f sec.") %"Network is constructed in" %lpr_bz.getCpuTime_Network())<<std::endl;
			std::cout<<(boost::format("   %-33s %8.2f sec.") %"Max Closure: Total CPU time" %total_MC_time)<<std::endl;
			std::cout<<(boost::format("   %-33s %8.2f sec.") %"Restricted LP: Total CPU time" %total_RLP_time)<<std::endl;
			
			std::cout<<(boost::format("%-36s %i") %"Number of iterations" %nIter)<<std::endl;
			std::cout<<(boost::format("%-36s %-28.6f") %"Objective value" %obj_val)<<std::endl;
			std::cout<<(boost::format("   %-33s %-9.6f %s") %"Gap" %((ub - obj_val)/obj_val) %"%")<<std::endl;
			std::cout<<(boost::format("\n%-36s %i") %"Number of periods" %nPeriods)<<std::endl<<std::endl;
			save_sol(argv[2], mcfactory.getOption(), status, obj_val, ub, X, Y, nPeriods, nIter, cpu_time, total_MC_time, lpr_bz.getCpuTime_Network(), total_RLP_time, 
											"_LPR_BZ.sol", "LP Solver: Bienstock-Zuckerberg Algorithm");
			
	}
	return 0;
}

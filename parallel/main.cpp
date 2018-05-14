
#include <iostream>
#include <mpi.h>
#include "MineProblem.h"
#include "../include/solver_interface.h"
#include "../include/daten.h"

//Allows disabling monitor thread for testing
#define MONITOR_THREAD

string input_file;
int diff_amt = 0;
int time_limit = 10; // a default amount of 10 minutes of run-time

// read the parameters provided as input
void read_parameters(int argc, char **argv){
  int iarg=1;
  while (iarg < argc)  {
    if (strcmp(argv[iarg],"-f")==0) {
      input_file = argv[++iarg];
    }
    else if (strcmp(argv[iarg],"-t")==0) { 
      time_limit=atoi(argv[++iarg]);
    }
    /*else if (strcmp(argv[iarg],"-maxiter")==0) {
      max_iter=atoi(argv[++iarg]);
    }
    else if (strcmp(argv[iarg],"-alg")==0){  
      alg_type=atoi(argv[++iarg]);
      }*/
    else if (strcmp(argv[iarg],"-diffamt")==0){  
      diff_amt=atoi(argv[++iarg]);
    }
    iarg++;
  }
}

int rank_alg_mapping(int rank){
  if(rank == 0) // ACO
    return 1;
  if(rank == 1) // BZ
    return 2;
  if(rank == 2) // LR - needs some debugging, not sure why it breaks during the second iteration
    return 1;
  if(rank == 3) // VLNS
    return 4;
  if(rank == 4) // VLNS fixed
    return 5;
  if(rank == 5) // Volume alg
    return 6;
  if(rank == 6) // LaPSO
    return 7;
  return 1;
}

int rank_alg_mapping_test(int rank){
  //  if(rank == 0) // ACO
  //  return 1;
  //if(rank == 1) // BZ
  //  return 4;
  //if(rank == 2) // LR - needs some debugging, not sure why it breaks during the second iteration
  //  return 4;
  //if(rank == 3) // VLNS
  //  return 4;
  //if(rank == 4) // VLNS fixed
  //  return 5;
  //if(rank == 5) // Volume alg
  //  return 6;
  //if(rank == 6) // LaPSO
  //  return 7;
  return 9;
  //return 3;
}

int main(int argc, char* argv[])
{
  if (argc < 2) 
  {
    std::cout << "Please provide a problem path argument\n";
    return 1;
  }

  MPI_Init(NULL, NULL);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  //Create a problem object
  //- reads the data, creates sub-problems, starts a worker thread
  //- sub-problems created per processor
  read_parameters(argc,argv);
  MineProblem problem(input_file.c_str(), rank, size);

  // Create a solver_interface instance here
  // Use it to invoke different solvers
  SolverInterface *si = new SolverInterface(*problem.data);
  // Set this for each run, just a constant for now
  //Get first problem
  problem.next();
  //double best_obj = -1e30;
  std::cout << "Time limit: " << time_limit << std::endl;

#ifdef MONITOR_THREAD
  //Start the worker task
  std::thread worker;
  worker = std::thread([&problem, &si, &rank]()
  {
    std::this_thread::sleep_for(std::chrono::milliseconds(1000*rank));
    std::cout << "Started solver thread: " << std::this_thread::get_id() << std::endl; 
#endif
    //Continue until all solved or timeout over (in minutes, TODO: pass limit as param)
    while (!problem.solved(true) && !problem.timeout(time_limit))
    //while (!problem.timeout(time_limit))
    {

      //Start solving active sub-problem
      // Specify some iterations to start with
      int iterations = 1000;
      // some dummy data as to test
      // Need the two things below to make the solver work (13/9/16)
      const int nB = problem.data->getNBlock();
      const int d_max = problem.data->getnDestination();
      //vector<int> X(nB,0);
      //vector<vector<double> > Y(nB, vector<double>(d_max,0.0));
      //problem.current->sol_int.x = X;
      //problem.current->sol_int.y = Y;
      problem.init_sol(nB, d_max);

      //int alg_type = rank+1; // this needs to change
      int alg_type = rank_alg_mapping_test(rank);
      //int alg_type = rank_alg_mapping(rank);
      std::cout << "Algorithm type: " << alg_type 
		<< ", rank: " << rank
		<< std::endl; 

      int shift = 2;
      int start = rank*shift; 
      std::cout << "Start: " << start 
		<< ", shift: " << shift
		<< std::endl; 
      si->callSolver(problem, diff_amt, iterations, alg_type, start, shift);
      if(alg_type == 9){
	Sol_Int *sol = si->getMipSol();
	Sol_Int * sol_master = new Sol_Int();
	problem.get_sol(*sol_master);
	delete sol_master;
      }
      //Solver returned with updated sub-problem, inform main thread/communicate with other procs...
      //if(problem.current->sol_int.obj > best_obj)
      //best_obj = problem.current->sol_int.obj;

      //Inform other processes we need to sync data
      problem.set_updated();
      //Sync data
      problem.sync();

      //Get next problem if current is solved
      //if (problem.solved())
      //  problem.next();
    }
#ifdef MONITOR_THREAD
    std::cout << "Solver thread done " << std::this_thread::get_id() << std::endl; 
  }
  );

  //Monitor on main thread
  problem.monitor();

  //Join solver thread
  worker.join();
#endif

  delete si;

  MPI_Finalize();
  return 0;
}

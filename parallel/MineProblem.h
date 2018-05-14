// Mine problem template:
//
//- Load problem data from path
//- Split into sub-problems
//- Each process (MPI)
//  - Starts a worker thread to run solver on own sub-problem set
//    (calls a black-box problem solver, multi-threaded via OpenMP or similar, no MPI code)
//  - Monitors worker thread in main thread, reporting back to master process
//- Master process
//  - Regularly queries status of all processes [idle/working]
//  - Report status to console
//  - As workers finish (enter idle state) retrieve their solutions
//  - Based on updated status:
//    - Stop all if solved state reached... 
//    - Stop workers whose problem tree is now redundant
//    - Redivide the problem and distribute new data
//    - Restart idle workers

//Questions:
// - How sub-problems creation should be managed and is this done on the root process or 
//   distributed, within worker class or before starting workers? Not sure how this will fit in
// - Solved state? A threshold presumably?
//Please edit as necessary with further methods or data required to encapsulate problem
#ifndef _MineProblem_
#define _MineProblem_

//Timeout in minutes
#ifndef TIMEOUT_LIMIT
#define TIMEOUT_LIMIT 1
#endif

#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <mpi.h>
#include <assert.h>

#include "../include/daten.h"
#include "../include/BranchNode.h"

//Container class for problem or sub-problem data
class ProblemData
{
public:
  BranchNode_info node;
  Sol_Int sol_int;
  Sol_Real sol_real;

  //Status
  bool solved;
  int assignedProc;

  ProblemData() : solved(false), assignedProc(-1)
  {
  }

    ProblemData(int nB,int nT,int nD,int nR) { init(nB,nT,nD,nR); }
    void init(int nB,int nT,int nD,int nR) {
	solved=false;
	node.init(nB,nT,nD,nR);
	sol_int.init(nB,nD);
	sol_real.init(nB,nT,nD);
    }
  ~ProblemData()
  {
  }

};

//Class to manage solving of a sub-problem
//- Worker thread will start and begin solving when run() is called
//- Monitor thread will sleep for a given time period, then request data from the worker thread
//- Status will then be returned
//- Main program will need to communicate data for all processes with MPI and determine 
//  whether to continue or update for each one
//- If update required, send the new problem data set
class MineProblem
{
public:
  int id; //Worker id (pass in mpi rank)
  int total; //Proc count

protected:
 public:  // delete this line - once BranchBound.cpp has been modified 
  int updated; //Updated flag
  std::atomic<bool> update_ready; //Updated received flag
  std::string problemPath;
  std::vector<ProblemData*> worklist;
  ProblemData* current;
  std::mutex mutex;
  std::chrono::system_clock::time_point start_time;

public:
  Daten* data;

  MineProblem(const char* path, int id, int count);
  ~MineProblem();

  bool timeout(int t_limit = 2)
  {
    auto current_time = std::chrono::high_resolution_clock::now();
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(current_time - start_time);

    if (minutes.count() >= t_limit)
    {
      std::cout << "Program has been running for " 
		<< minutes.count() << " minutes, aborting" 
		<< ", Limit: " << t_limit << " minutes." << std::endl;
      return true;
    }
    return false;
  }

  //Problem sub-division data/status
  int openSubProblems() { return worklist.size(); }
  int solvedSubProblems()
  {
    std::lock_guard<std::mutex> guard(mutex);
    //Iterate worklist and count processed
    int count = 0;
    for (ProblemData* problem : worklist)
      if (problem->solved) count++;
    return count;
  } 

  bool solved(bool all=false) 
  {
    std::lock_guard<std::mutex> guard(mutex);
    if (!current) return true; //No active problem

    //Checking if current problem solved
    if (!all)
      return current->solved;

    //Checking if all problems solved
    for (ProblemData* problem : worklist)
      if (!problem->solved) return false;
    return true;
  }

  void setSolved()
  {
    std::lock_guard<std::mutex> guard(mutex);
    if (current) current->solved = true; //Set active problem to solved
  }

  void split();
  void addWork(ProblemData *child)  { child->solved = false; worklist.push_back(child); } 

  //Get next problem to solve
  void next();

  //Monitor solver
  void monitor();

  void setUpdated();
  bool sync();


    // TODO: we need to be able to get the branch node & modify it for
    //       preprocessing. at the moment this works ok as is but thread
    //       issues might arise (?)    
  BranchNode_info& get_BranchNode(bool lock=false)
  {
    //Return a reference for now, once we start modifying/updating this, need to:
    // - make a copy first, or
    // - pass lock=true and call release_BranchNode() when done updating
    if (lock) mutex.lock();
    //Return a reference so can be modified, when done modifying unlock mutex
    return current->node;
  }

  void release_BranchNode()
  {
    //Unlock the mutex to allow update
    mutex.unlock();
  }

  void set_updated()
  {
    std::lock_guard<std::mutex> guard(mutex);
    //std::cout << id << " Finished step, OBJ = " << current->sol_int.obj << std::endl;
    updated = true;
  }

  bool get_updated()
  {
    std::lock_guard<std::mutex> guard(mutex);
    return updated > 0;
  }

  void init_sol(int nB, int d_max)
  {
    std::lock_guard<std::mutex> guard(mutex);
    current->sol_int.init(nB, d_max);
  }

  void update_sol(const Sol_Real &sol)
  {
    if(current->sol_real.obj > sol.obj) { // better upper bound
      std::lock_guard<std::mutex> guard(mutex);
      current->sol_real = sol;
    }
  }
  void update_sol(const Sol_Int &sol)
  {
    if(current->sol_int.obj < sol.obj) { // better lower bound
      std::lock_guard<std::mutex> guard(mutex);
      current->sol_int = sol;
      //std::cout << "\t\tNew best found: updating problem .." << std::endl;
    }
  }


  void update_sol(const int &nT, const double &obj, const std::vector<int> &x, const std::vector< std::vector<double> > &y)
  {
    std::lock_guard<std::mutex> guard(mutex);
    current->sol_int.obj = obj;
    current->sol_int.nT = nT;
    current->sol_int.x = x;
    current->sol_int.y = y;
  }

  void get_sol(int &nT, double &obj, std::vector<int> &x, std::vector< std::vector<double> > &y)
  {
    std::lock_guard<std::mutex> guard(mutex);
    obj = current->sol_int.obj;
    nT = current->sol_int.nT;
    x = current->sol_int.x;
    y = current->sol_int.y;
  }

  void get_sol(Sol_Int &sol)
  {    std::lock_guard<std::mutex> guard(mutex);
       sol = current->sol_int;
  }
  void get_sol(Sol_Real &sol)
  {    std::lock_guard<std::mutex> guard(mutex);
       sol = current->sol_real;
  }
  double get_ub() {
      std::lock_guard<std::mutex> guard(mutex);
      return current->sol_real.obj;
  }
  double get_obj()
  {
    std::lock_guard<std::mutex> guard(mutex);
    return current->sol_int.obj;
  }
};

#if 0
//Now using solver_interface.h version...

//Solver interface, extend this to implement a particular solver algorithm
class Solver
{
 public:
  ProblemData* solving;

  Solver() {}

  virtual void solve(MineProblem& problem) = 0;
};

class TestSolver : public Solver
{
 public:

  TestSolver() : Solver() 
  {
    //Init RNG with MPI rank so each proc gets different results
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    srand(rank); //Use MPI rank here to trigger different search per proc
  }

  virtual void solve(MineProblem& problem)
  {
    solving = problem.current;
    solving->node.ub = 1000000;
    //Loop until solved or interrupted
    while (1)
    {
      //Solver implementation must lock this mutex before using "problem" data
      //however it must not remain locked to give the monitor thread a chance to 
      //check for updates, so here we are locking it at the start of each solver iteration
      problem.mutex.lock();

      //!!!!Current issue with this:
      //Each proc calling update regularly holds up all processes until they make the call
      //Needs to be called frequently or use some other method
      //
      //Updated by another process?
      problem.mutex.unlock();
      //Call the update function, all procs must call this to get updated data
      //  After the update we may want to break and abandon this sub-problem...
      //  currently signalled by returning false, or could decide here within solver
      if (!problem.update())
      {
        //Aborting, flag as solved so skipped in next()
        solving->solved = true;
        break;
      }
      problem.mutex.lock();
      
      double value = ((double)rand()) / (double)RAND_MAX;
      //Increasing obj (profit)
      solving->sol_int.obj += value;

      //Lets just say for now it's solved if obj greater than 100.0
      //Solved?
      if (solving->sol_int.obj > 100.0)
      {
        solving->solved = true;
        break;
      }
      problem.mutex.unlock();
      std::cout << "OBJ " << solving->sol_int.obj << std::endl;
    }
    std::cout << "Exiting solver...\n";

  }
};

#endif

#endif

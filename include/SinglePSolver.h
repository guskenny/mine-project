#ifndef __SinglePSolver_H__
#define __SinglePSolver_H__

#include "SinglePModel.h"
#include "MaxClosureFactory.h"
#include "BranchNode.h"
#include "Preprocess.h"
#include <boost/format.hpp>
#include "QOL/QolMIP.h"
#include "QOL/QolColFormulation.h"
#include "QOL/CplexFormulation.h"
#include "QOL/GurobiFormulation.h"
//#include "QOL/CPBeamACO.h" // needs compilation with boost::thread
#include "QOL/CpuTimer.h"
#include <my_math.h>
#include <daten.h>
#include <SinglePModel.h>
#include <cmath>
#include <map>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <algorithm>
#include "SettingsHandler.h"
#include "SolutionMerger.h"
#include "RandomSearch.h"
#include "LocalSearch.h"
#include "MergeSolver.h"
#include "MergeSolverSimple.h"
#include "MergeSolverCompact.h"
#include "ConeMiner.h"
//#include <math.h>
// #include <ncurses.h>
#include <sstream>
#include <queue>

typedef std::pair<double, int> p_mined_pair;

#define CPLEX_T 1
#define GUROBI_T 2

class SinglePSolver{
  private:
    SinglePModel *probModel;
    MaxClosure_Base *MCSolver;
    std::vector<double> residualProfit;
    char MCType;
    int solverType;
    bool solveRelaxed;
    SettingsHandler sh;

  public:
    SinglePSolver(const Daten &prob);
    SinglePSolver(int argc,const char **argv);

    void initModel(qol::MIPSolver &mip,
        std::vector<qol::Variable> &x,
        //std::vector<std::vector<qol::Variable> > &y,
        std::vector<Block> * blocks,
        std::vector<double> &profitModifier);

    void computeResUse(Sol_Int &sol);

    int solve();

    void getSeeds(BranchNode_info blank_info,std::vector<Sol_Int> &seeds);
    void computeUPIT(BranchNode_info &base_info, std::vector<int> &include);

    void saveSols(const std::vector<Sol_Int> &sols, std::string path_name);
    bool loadSols(std::vector<Sol_Int> &sols);

    int forkMergeSolve();
    int serialMergeSolve();

    double singlePSolve(const BranchNode_info &init_sol, Sol_Int &sol);

    void doMerge(SolutionMerger &sm, const std::vector<Sol_Int>&sols, const std::vector<int> &include, Sol_Int &init_sol, Sol_Int &merged_sol, std::ofstream &red_data);

    void setSolverType(int type) {solverType = type;};

    ~SinglePSolver(){};

};

#endif

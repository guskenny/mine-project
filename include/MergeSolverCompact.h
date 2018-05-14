#ifndef __MergeSolverCompact_H__
#define __MergeSolverCompact_H__

#include "Preprocess.h"
#include <boost/format.hpp>
#include "QOL/QolMIP.h"
#include "QOL/QolColFormulation.h"
#include "QOL/CplexFormulation.h"
#include "QOL/GurobiFormulation.h"
#include "QOL/CpuTimer.h"
#include <daten.h>
#include <math.h>
#include <map>
#include <vector>
#include <set>
#include "graph.h"
#include "SettingsHandler.h"
#include "BranchNode.h"
#include "QOL/CpuTimer.h"
#include "SinglePModel.h"

class MergeSolverCompact{
  private:
    SinglePModel *probModel;
    SettingsHandler sh;
    Sol_Int best_sol;
    int nB;
    int nG;
    int nR;
    int nT;
    double rate;
    std::vector<std::vector<int> > x_map;
    std::vector<qol::Variable> x;
    std::vector<std::vector<int> > groups;
    std::vector<int> group_map;
    std::vector<int> fixed;
    Graph merge_graph;
    std::vector<int> include;
    std::ofstream *red_data;

  public:
    bool rel_gap;
    const int BACKWARD = 0;
    const int FORWARD = 1;
    MergeSolverCompact(const SettingsHandler sh, SinglePModel *base_model, Sol_Int &best_sol, std::vector<std::vector<int> > &groups, std::vector<int> &group_map, std::vector<int> &fixed, const std::vector<int> &include,std::ofstream *red_data);

    void solve(Sol_Int &sol);

    void initMergeGroupModel(qol::MIPSolver &mip);
    void initMergeSimpleModel(qol::MIPSolver &mip);

    void extractSimpleSol(qol::MIPSolver &mip, Sol_Int &sol);
    void extractGroupSol(qol::MIPSolver &mip, Sol_Int &sol);

    ~MergeSolverCompact(){};

};

#endif

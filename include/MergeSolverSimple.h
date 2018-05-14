#ifndef __MergeSolverSimple_H__
#define __MergeSolverSimple_H__

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

class MergeSolverSimple{
  private:
    SinglePModel *probModel;
    SettingsHandler sh;
    Sol_Int best_sol;
    int nB;
    int nG;
    int nR;
    int nT;
    double rate;
    std::vector<std::vector<qol::Variable> > x;
    std::vector<std::vector<qol::Variable> > y;
    std::vector<std::vector<int> > groups;
    std::vector<int> group_map;
    std::vector<int> fixed;
    Graph merge_graph;

  public:
    bool rel_gap;
    const int BACKWARD = 0;
    const int FORWARD = 1;
    MergeSolverSimple(const SettingsHandler sh, SinglePModel *base_model, Sol_Int &best_sol, std::vector<std::vector<int> > &groups, std::vector<int> &group_map, std::vector<int> &fixed);

    void solve(Sol_Int &sol);
    void initMergeModel(qol::MIPSolver &mip);
    void initMergeSimpleModel(qol::MIPSolver &mip);

    ~MergeSolverSimple(){};

};

#endif

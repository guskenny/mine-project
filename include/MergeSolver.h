#ifndef __MergeSolver_H__
#define __MergeSolver_H__

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

class MergeSolver{
  private:
    SettingsHandler sh;
    BranchNode_info merged;
    int nB;
    int nG;
    int nR;
    int nT;
    double rate;
    std::vector<std::vector<int> > groups;
    std::vector<int> group_map;
    std::vector<double> group_profits;
    std::vector<std::vector<double> >group_resources;
    Graph merge_graph;

  public:
    const int BACKWARD = 0;
    const int FORWARD = 1;
    MergeSolver(const SettingsHandler sh, SinglePModel *base_model, BranchNode_info &merged, std::vector<std::vector<int> > &groups, std::vector<int> &group_map);

    void initMergeModel(SinglePModel *base_model);

    ~MergeSolver(){};

};

#endif

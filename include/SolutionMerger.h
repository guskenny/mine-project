#ifndef __SolutionMerger_H__
#define __SolutionMerger_H__

#include <daten.h>
#include <math.h>
#include <map>
#include <vector>
#include <set>
#include <queue>
#include "graph.h"
#include "SettingsHandler.h"
#include "BranchNode.h"
#include "QOL/CpuTimer.h"
#include "SinglePModel.h"

class SolutionMerger{
  private:
    SettingsHandler sh;

  public:
    SolutionMerger(const SettingsHandler sh) : sh(sh) {std::cout << "solution merger initialised" << std::endl;};

    void fullMerge(const std::vector<Sol_Int>&sols,const std::vector<int> &include, std::vector<int> &fixed, std::vector<std::vector<int> > &groups, std::vector<int> &group_map);

    void simpleMerge(const std::vector<Sol_Int>&sols, const std::vector<int> &include, std::vector<int> &fixed);
    void mergeCPIT(const std::vector<Sol_Int>&sols, std::vector<int> &fixed, std::vector<std::vector<int> > &groups, std::vector<int> &group_map);
    void mergePCPSP(const std::vector<Sol_Int> &sols, BranchNode_info &merged);
    void writePCPSPmerged2D(const std::vector<Sol_Int> &sols);
    void writePCPSPmerged3D(const std::vector<Sol_Int> &sols);

    ~SolutionMerger(){};

};

#endif

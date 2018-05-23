#ifndef __LocalSearch_H__
#define __LocalSearch_H__

#include <daten.h>
#include <math.h>
#include <random>
#include <utility>
#include <fstream>
#include "SetObj.h"
#include "graph.h"
#include "SinglePModel.h"
#include "SettingsHandler.h"
#include "BranchNode.h"
#include <boost/format.hpp>
#include "QOL/QolMIP.h"
#include "QOL/QolColFormulation.h"
#include "QOL/CplexFormulation.h"
#include "QOL/GurobiFormulation.h"
//#include "QOL/CPBeamACO.h" // needs compilation with boost::thread
#include "QOL/CpuTimer.h"

typedef std::pair<double, int> p_mined_pair;

class LocalSearch{
  private:
    // const int BACKWARD = 0;
    // const int FORWARD = 1;
    // const int NUM_DIRS = 2;

    SettingsHandler sh;
    std::mt19937 rng;
    SinglePModel *model;
    std::vector<int> include;

  public:
    LocalSearch(const SettingsHandler sh, SinglePModel *model, const std::vector<int> &include) : sh(sh), model(model), include(include) {
      std::cout << "LocalSearch initialised" << std::endl;
      std::random_device r;
      std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
      rng = std::mt19937(seed);
    };

    const int getDirection();

    bool computeCone(const Sol_Int &sol, const int &cone_tip,
      SetObj &cone, long double &cone_profit, std::vector<long double> &cone_res_use,
      const int &direction, int &cone_depth);

    void swapConePeriod(Sol_Int &sol, const int period, const int dir,
      const SetObj &cone, const long double cone_profit, const std::vector<long double> &cone_res_use, std::vector<SetObj> &period_blocks);

    void SAnoPeriod(Sol_Int &sol);

    void repairSolution(Sol_Int &sol);

    void swapWalk(Sol_Int &sol);
    void goodSwap(Sol_Int &sol);

    void runTests(Sol_Int &sol);

    void runSingleTests(Sol_Int &sol);

    void displayGraph(int width, Sol_Int &sol,
      std::vector<SetObj> &bounds, const SetObj &cone);
    void genTestGraph(int nB, int width, Graph &graph);

    ~LocalSearch(){};

};

#endif

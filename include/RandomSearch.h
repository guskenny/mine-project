#ifndef __RandomSearch_H__
#define __RandomSearch_H__

#include <daten.h>
#include <math.h>
#include <random>
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

class RandomSearch{
  private:
    SettingsHandler sh;
    std::mt19937 rng;

  public:
    RandomSearch(const SettingsHandler sh) : sh(sh) {
      std::cout << "RandomSearch initialised" << std::endl;
      std::random_device r;
      std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
      rng = std::mt19937(seed);
    };
    double randomSearch(const BranchNode_info &probInfo, Sol_Int &sol, SinglePModel *model);
    double betterRandomSearch(const BranchNode_info &probInfo, Sol_Int &sol,SinglePModel *model);

    ~RandomSearch(){};

};

#endif

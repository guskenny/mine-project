#ifndef __GRASPSolver_H__
#define __GRASPSolver_H__

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

class GRASPSolver{
  private:
    // const int BACKWARD = 0;
    // const int FORWARD = 1;
    // const int NUM_DIRS = 2;
    int nB;
    int nR;
    int nT;
    double rate;
    SettingsHandler sh;
    SinglePModel *probModel;
    std::vector<int> include;

  public:
    GRASPSolver(const SettingsHandler sh, SinglePModel *model, const std::vector<int> &include) : sh(sh), probModel(model), include(include) {
      std::cout << "GRASPSolver initialised" << std::endl;
      nB = model->getNBlock();
      nT = model->getNPeriod();
      nR = model->getnResources();
      rate = model->getDiscountRate();
    };

    void solveWindow(Sol_Int &sol, int t_0);
    void initWindowModel(qol::MIPSolver &mip, std::vector<std::vector<qol::Variable> > &x, const std::vector<int> &in_window, const std::vector<int> &window_map,Sol_Int &sol, int t_0,const std::vector<int> &fixed);
    
    ~GRASPSolver(){};

};

#endif

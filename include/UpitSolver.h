#ifndef UpitSolver_H
#define UpitSolver_H

#include "daten.h"
#include "MaxClosure_BoostMaxFlow_BK.h"

class UpitSolver {
  private:
    size_t nb;
    MaxClosure_Base * upit;
    
  public:  
    UpitSolver(const Daten &prob);
    ~UpitSolver(){};
    int solve();
    int getFixed(std::map<int,int> &fixed);
    int getClosure(std::vector<int> &closure);
    void getResidualProfit(std::vector<double> &residualProfit);
};

#endif

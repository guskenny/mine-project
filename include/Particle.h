#include "QOL/lagrangian/LaPSO.hpp"
#include "QOL/lagrangian/Random.h"
#include "QOL/CpuTimer.h"
#include "CumulativeModel.h"
#include "MaxClosure_Base.h"

#define SCALE 1e6 // divide constraints and profits by SCALE

class PCPSPparticle : public LaPSO::Particle {
public:
    //PCPSPparticle() : mc(0){}
    PCPSPparticle(int nVertex,int nResources,MaxClosure_Base *m)
	: LaPSO::Particle(nVertex,nResources),mc(m) {}
    ~PCPSPparticle() { if(mc!=0) delete mc; }
    MaxClosure_Base *mc;		// max closure solver - need one for each thread/particle
    LaPSO::Uniform rand;		// independent random number generator
    //double lagConst;			// current lagrangian constant
};

// for each vertex a list of resource and amount consumed vertex is included
// (ie the inverse of CumulativeModel.res
class ResourceUse : std::vector<std::vector<std::pair<int,double> > >{
public:
    ResourceUse(const CumulativeModel &prob) ;
    // can the vertex v be included in the solution given residual amount of
    // resource left
    bool canFit(CumulativeModel::Vertex v,const std::vector<double> &residual) const;
    // update residual by subtracting out resource usage when vertex v is included
    void update(CumulativeModel::Vertex v,std::vector<double> &residual) const;
    // undoUpdate - negate the effects of an update
    void undoUpdate(CumulativeModel::Vertex v,std::vector<double> &residual) const;
};


class PCPSPpsoHooks : public LaPSO::UserHooks {
public:
    CumulativeModel &prob;
    ResourceUse resUse;
    double bestCost;
    LaPSO::IntVec bestSoln;
    qol::CpuTimer timer;
    int nMCsolves;
    PCPSPpsoHooks(CumulativeModel &model) :
	prob(model),resUse(model),bestCost(-1e99), bestSoln(model.graph.getNumNodes()),
	nMCsolves(0) {}
    LaPSO::Status reducedCost(const LaPSO::Particle &p_, LaPSO::DblVec &redCost);
    // LaPSO solves max_dual min_x -profit * x - dual (resLim - resUse(x))
    // hence redCost = -1 * prob.calcRedCost()
    LaPSO::Status solveSubproblem(LaPSO::Particle &p_);
    LaPSO::Status fixConstraint(const int res,const LaPSO::Particle &p_,LaPSO::SparseVec &feas);
    LaPSO::Status heuristics(LaPSO::Particle &p_) ;
    LaPSO::Status updateBest(LaPSO::Particle &p_);
};

// main function for Lagrangian particle swarm or Volume algorithm
// runs LaPSO if -L argument is included else with -V runs volume algorithm
// returned problem should be de-allocated when all done
LaPSO::Problem *runParticleMain(CumulativeModel &prob,int argc,const char **argv); 

    

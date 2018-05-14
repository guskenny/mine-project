/* LP Relaxation: Open Pit Mining(PCPSP)
 * author: Davaa Baatar
 * date: 19.06.2016
*/
#ifndef __RLP_H__
#define __RLP_H__
#include <util.h>
#include "BranchNode.h"
#include "MineProblem.h"
#include "CumulativeModel.h"

/// Combining solutions with this module:
///  create Partition and setFixed() with splitAll=false
/// for each (heuristic/lagrangian/...) solution call split()
/// create ResLP(partition) and call solveExact()

// Partition[i] is set of vertices in cummulative graph that make up the set
class Partition : public std::vector<std::vector<int> > {
protected:
    int n;		    // number of vertices
    std::vector<int> fixed; // for each set whether it is fixed to zero or one (-1=free)
public:
    enum FixedValue { Free=-1, FixedZero=0,FixedOne=1};
    // set up parititon with n elements
    // if provided the
    Partition(const CumulativeModel &model,const BranchNode_info *node=0) :
	std::vector<std::vector<int> >(1),n(model.graph.getNumNodes()), fixed(1,-1)
	{ if(node) setFixed(model,*node); }
    Partition(int _n) : std::vector<std::vector<int> >(1),n(_n), fixed(1,-1)	
	{ back().resize(n); for(int i=0;i<n;++i) back()[i] = i; }
    // setting fixed resets partition into (up to) 3 sets, fixed zero/one/free
    // if splitAll=true it creates a partition with no aggregation that isn't
    // required by the BranchNode_info.
    void setFixed(const Daten &model,const BranchNode_info &node,bool splitAll=false);
    bool isFixed(int s) const {return fixed[s]!=-1;}	// is partition fixed?
    int fixedVal(int s) const {return fixed[s];}  // all elements = fixed val.
    void fix(int s,int val) { fixed[s]=val;}	  // fix partition to value
    // change partition so that elements in a set have the same value for cut
    // returns the number of sets split (= increase in size of Partition
    int split(const std::vector<int> &cut); // cut must be length n.    
    // update partitition in a way that allows this solution to be reproduced
    // can be done multiple times 
    int split(const CumulativeModel &model,const Sol_Int &soln); 
};

class ResLP{
public:
    // Cpart[setNo] -> vector or vertices in the cummulative graph that make up the set (all the same value) (INPUT)
    ResLP(CumulativeModel & data,const Partition& Cpart);
    ~ResLP();
    // after solve, solution stored in references provided
    // pi[time][resource]  - dual value (output)
    // sol[vertex] -> fraction  (output)
    void solve(std::vector<double> & sol,std::vector<std::vector<double> > & pi);
    // use mip solver to get exact solution - post any integer solutions via update
    // in future best may also be used as input to give CPLEX/Gurobi a
    // starting solution
    void solveExact(MineProblem &prob,Sol_Int &best,double timeLimit=3600);
    void solve_LargeLP(std::vector<double> & sol,std::vector<std::vector<double> > & pi);  // unused, for testing only
    double getObjVal();  // LP value
    double getTime();   // CPU time
    bool getStatus();  // true = solved, false = infeasible/failed

private:
    CumulativeModel  & _data;
    //std::vector<double> & _sol;
    double   _objVal;
    bool _status;
    //std::vector<std::vector<double> > & _pi;
    const Partition & _Cpart;
    double _t_cpu;
};
#endif

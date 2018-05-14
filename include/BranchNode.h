// type definitions
#ifndef __BRANCH_NODE_H__
#define __BRANCH_NODE_H__

#include <array>

// test

//-------------------------------------  Tree  --------------------------------------------
struct BranchNode_info{

    BranchNode_info() : ub(1e99) {}
    // initalise with number of blocks,time periods,destinations, and resources
    BranchNode_info(int nB,int nT,int nD,int nR) { init(nB,nT,nD,nR); }
    void init(int nB,int nT,int nD,int nR) {
	time.resize(nB);
	dest.resize(nB);
  fixed.resize(nB);
	mine_block.resize(nB);
	ub=1e99;
	w_start.resize(nT);
	for(int b=0;b<nB;++b){
	    time[b][0]=0;
	    time[b][1]=nT;
      fixed[b][0]=-1;
      fixed[b][1]=nT;
	    dest[b].resize(2,0.0);
	    mine_block[b] = true; // available to mine
	}
	for(int t=0;t<nT;++t)
	    w_start[t].resize(nR,0.0);
    }

    int fixedCnt() const {
	const int nT = w_start.size();
	int cnt=0;
	for(size_t b=0;b<time.size();++b) cnt += time[b][0]+(nT-time[b][1]);
	return cnt;
    }

  // Block b can only be mined during period time[b][0],...,time[b][1]-1
  // if time[b][0]=time[b][1] it can't be mined at all
  // if time[b][0]=time[b][1]-1 it must be mined in this period

  //TODO: add extra branchnode data that takes fixed time periods because it is not the same as restricted time interals
  std::vector<std::array<int, 2> > fixed;
  std::vector<std::array<int, 2> > time; // for each block a restricted time interval
  std::vector<std::vector<double> > dest; // amount to send to each destination (can be ignored by solver)
  // the following information may be passed from parent nodes during a branch and bound process:
  double ub; // relaxed bound (initial best bound)
  //w_start[t][r] is the dual (lagrange multiplier) for period t, resource r (optional)
  std::vector< std::vector<double> > w_start;     // warm start(Lagrangian multipliers)
  // Whether a block must be mined or not
  // If yes, use time (above) to determine what period.
  // if no, mining during time interval is optional
  std::vector<bool> mine_block;  // mine_block[b] -> include b in ultimate pit
};

//-------------------------------------  Solution  --------------------------------------------
struct Sol_Int{ // feasible (integral) solution
    Sol_Int() : obj(0.0), nT(0), nR(0) {}
    Sol_Int(int nB, int nT) { init(nB,nT); }
    Sol_Int(int nB, int nD, int r_max, int t_max) {init(nB, nD, r_max, t_max); }
    void init(int nB, int _nT) {
        obj = 0.0;
        nT = _nT;
        nR = 2;
        res_use.resize(_nT, std::vector<long double>(2, 0.0));
        x.resize(nB, _nT);
        y.resize(nB, std::vector<double> (1, 0.0));
    }
    void init(int nB, int nD, int r_max, int t_max){
      init(nB, t_max);
      res_use.resize(t_max, std::vector<long double>(r_max, 0.0));
      nR = r_max;
      nT = t_max;
    }
    // Sol_Int & operator=(const Sol_Int & rhs) // not needed as C++ defines by default

    std::vector<int> x;                    // block => period; nT if not excavated
    std::vector< std::vector<double> > y;  // y[b][d] = fraction of block b sent to the destinations d;
    std::vector<std::vector<long double> > res_use;
    int nT;                                // number of periods
    int nR;
    double obj;                            // obj value
    int numMined() const {
	int k=0; for(int t : x){ if(t<nT) ++k;} return k;
    }
};

// Note: Sol_Real is _VERY_ large due to y matrix (blocks * periods * destinations ~ millions of doubles)
//       Expect not to be passing this between compute nodes.
struct Sol_Real{
    Sol_Real() : obj(1e99) {}
    Sol_Real(int nB, int T, int nD) { init(nB, T, nD); }
    Sol_Real(std::vector<std::vector<double> > & Xval, std::vector< std::vector< std::vector<double> > > & Yval, int p, double objval){
        x = Xval;
        y = Yval;
        nT = p;
        obj = objval;
    }
    void init(int nB, int T, int nD) {
        x.resize(nB, std::vector<double> (T, 0.0));
        y.resize(nB, std::vector< std::vector<double> > (T, std::vector<double> (nD, 0.0)));
        obj = 1e99;		/* infinite upper bound */
        nT = 0;
    }
    std::vector<std::vector<double> > x;                     // block -> time -> value
    std::vector< std::vector< std::vector<double> > > y;     // block->time-> destination => value
    int nT;                                                   // number of periods
    double obj;                                              // obj value
};
#endif

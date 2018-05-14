/* LP Relaxation: Open Pit Mining(PCPSP)
 * Bienstock-Zuckerberg algorithm(Initial Version)
 * author: Davaa Baatar
 * date: 19.06.2016
*/
#ifndef __BZ_H__
#define __BZ_H__
#include "MaxClosure_Base.h"
#include "CumulativeModel.h"
#include "MaxClosureFactory.h"
#include "../QOL/CpuTimer.h"
#include <boost/format.hpp>
#include "rlp.h"
#include "BranchNode.h"
#include "MineProblem.h"

class BZ{
	public:
                BZ( CumulativeModel &prob, MaxClosureFactory mcfactory, double timeLimit = -1,
		    const BranchNode_info *branch=0); // if provided is copied into class
		~BZ();
		//method
                bool solve(MineProblem *best=0); // update best solution if best!=0
		
		// enquiry
		double getCpuTime();
		double getCpuTime_MC();
		double getCpuTime_RLP();
		double getCpuTime_Network(); // network construction time
		int get_nIter();
		int getStatus();
		
		double getObjValue();
		double getUB();
		std::vector<double> & getSolution();
		BranchNode_info &getBranchNode() {return _branch;}
		const BranchNode_info &getBranchNode() const { return _branch; }
		
	private:
		// members
		qol::CpuTimer _timer;
		double _timelimit;
		int _nIter;
		double _cpu_time;
		double _cpu_mc;
		double _cpu_rlp;
		double _cpu_network;
		int _status;
		
		CumulativeModel & _prob;
		BranchNode_info _branch;
		MaxClosure_Base * _mc;
		
		double _obj_const;
		double _obj_val;	// lower bound from restricted LP
		double _ub;			// UB from Lagrangian relaxation
		
		int _nVertices;
		std::vector<double> _sol_z;
		std::vector<int> _maxCut;	// 0/1 
		
		Partition _Cpart; // the partitioning: partition index -> the set of vertices
		std::vector<double> _lambda; 	
		std::vector<double> _supply;  // NOTE: first element is for the dummy node
		
		// methods
		int updatePartitioning(); // return 1 if stop criteria holds
		void reindex_Multipliers(std::vector<std::vector<double> > & lambda);
		void restore_Multipliers(const std::vector<std::vector<double> > & lambda);

};
#endif

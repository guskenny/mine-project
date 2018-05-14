// type definitions
#ifndef __TYPEDEF_H__
#define __TYPEDEF_H__
//-------------------------------------  Tree  --------------------------------------------
struct BranchNode_info{
	BranchNode_info(int i, int period, bool btype, std::vector< std::vector<double> > v, double ubound){
		_id = i;
		_t = period;
		_type = btype;
	}
	BranchNode_info & operator=(const BranchNode_info & rhs){
		_id = rhs._id;              							// block id
		_t = rhs._t;			     							// period
		_type = rhs._type;		 								// branch types: type 0 = at fixed time,	type 1 = defined w.r.t. "earlier than" including t
		return *this;
	}
	
	int _id;              							// block id
	int _t;			     							// period
	bool _type;		 								// branch types: type 0 = at fixed time,	type 1 = defined w.r.t. "earlier than" including t
};

struct Branch_history{
	Branch_history(){}
	Branch_history(std::vector< BranchNode_info > & ones, std::vector< BranchNode_info > & zeros){
		_blocks_one = ones;
		_blocks_zero = zeros;
	}
	Branch_history & operator=(const Branch_history & rhs){
		_blocks_one = rhs._blocks_one;
		_blocks_zero = rhs._blocks_zero;
		return *this;
	}
	std::vector< BranchNode_info > _blocks_one;	
	std::vector< BranchNode_info > _blocks_zero;	
};

struct Branch_Node{
	int _id;
	int _depth;
	double _ub;
	BranchNode_info _info;
	Branch_history _history;
	std::vector< std::vector<double> > _w_start; 	// warm start(Lagrangian multipliers)
};
//-------------------------------------  Solution  --------------------------------------------
struct Sol_Int{
	Sol_Int(){}
	Sol_Int(int nB, int nD){
		X = std::vector<int> (nB, -1);
		Y = std::vector< std::vector<double> > (nB, std::vector<double> (nD, 0.0));
		obj = 0.0;
		t = 0;
	}
	Sol_Int & operator=(const Sol_Int & rhs){
		X = rhs.X;
		Y = rhs.Y;
		t = rhs.t;
		obj = rhs.obj;
		return *this;
	}
	std::vector<int> X;					   // block => period; -1 if not excavated
	std::vector< std::vector<double> > Y;  // amount sent to the destinations; block -> destination => amount
	int t;								   // number of periods
	double obj;							   // obj value
};

struct Sol_Real{
	Sol_Real(){}
	Sol_Real(int nB, int T, int nD){
		X = std::vector<std::vector<double> > (nB, std::vector<double> (T, 0.0));
		Y = std::vector< std::vector< std::vector<double> > > (nB, std::vector< std::vector<double> > (T, std::vector<double> (nD, 0.0))); 
		obj = 0.0;
		t = 0;
	}
	Sol_Real(std::vector<std::vector<double> > & Xval, std::vector< std::vector< std::vector<double> > > & Yval, int p, double objval){
		X = Xval;
		Y = Yval;
		t = p;
		obj = objval;
	}
	Sol_Real & operator=(const Sol_Real & rhs){
		X = rhs.X;
		Y = rhs.Y;
		t = rhs.t;
		obj = rhs.obj;
		return *this;
	}
	std::vector<std::vector<double> > X; 					// block -> time -> value
	std::vector< std::vector< std::vector<double> > > Y; 	// block->time-> destination => value
	int t;								   					// number of periods
	double obj;							   					// obj value
};
#endif

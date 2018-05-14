/* LP Relaxation: Open Pit Mining(PCPSP)
 * author: Davaa Baatar
 * date: 19.06.2016
*/
#ifndef __MAXCLOSURE_H__
#define __MAXCLOSURE_H__
#include <daten.h>
#include <string>
//#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h> 
#include <cstdlib>
#include <math.h> 
#include <time.h>
#include <util.h>
class MaxClosure{
public:
	MaxClosure(Daten & data);
	~MaxClosure();
	void solve(double * demand);
	double getObjVal();
	double getTime();
	bool updatePartitioning(std::vector<std::vector<int> > & Cpart);
private:
	Daten & _data;
	int * _tails;
	int * _heads;
	double * _obj;
	std::vector<int> _sol; // index -> {0,1}
	double _objVal;
	int _nArcs;
	int _nVertices;
	CPXENVptr _env;
	CPXNETptr _net;	
	int _solstat;
	double _t_cpu;
};
#endif

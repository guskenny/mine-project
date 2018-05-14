/* LP Relaxation: Open Pit Mining(PCPSP)
 * author: Davaa Baatar
 * date: 19.06.2016
*/
#ifndef __UTUL_CPP__
#define __UTUL_CPP__
#include <daten.h>
#include <vector>
#include <string>
#include <boost/format.hpp>
#include "BranchNode.h"
const double eps_tol = 1e-6;
const double res_eps_tol = 1e-6;
const double obj_eps_tol = 1e-2; 
// const double eps_tol = 1;
//#include <iomanip>
void info_data(Daten data, char * filename);
void usage(char ** argv );
void save_sol(char * filename, char opt, bool status, const double obj_val,
	const double ub, const std::vector<std::vector<double> > & X,
	const std::vector< std::vector< std::vector<double> > > & Y, const int nPeriods,
	const int nIter, const double total_time, const double total_MC_time, const double cpu_network,
	const double total_RLP_time, std::string attr, std::string text);
void save_sol_int(char * filename, bool status, const Sol_Int & sol, const double obj_val, const double total_time, std::string attr, std::string text);



/* evaluate objective coefficients w.r.t. lambda */
double evaluate_const(Daten & data, std::vector<std::vector<double> > & lambda);
void evaluate_supply(Daten & data, std::vector<std::vector<double> > & lambda, double * supply);

/* mappings: (btd) -> index;  index->(btd) */
int convert_triplet_to_index(const Daten & data, int b, int t, int d);
void convert_index_to_triplet(const Daten & data, int index, int & b, int & t, int & d);

/* get opt. solution of the PCPSP(x,y) from opt. sol. of PCPSP(z) */
void convert_Z_to_XY(Daten & data, std::vector<double> & sol_z, std::vector<std::vector<double> > & X, std::vector< std::vector< std::vector<double> > > & Y, int &nPeriods);

// verify solution return 1 if failed, 0 otherwise
bool verify(const Daten & data, const Sol_Int & sol);
bool verify(const Daten & data, const Sol_Real & sol);

/* evaluate objective coefficients of the Restricted LP */

/* evaluate resource constraint coefficients of the Restricted LP */

/* evaluate the objective value */
double evaluateObjValue(const Daten & data, Sol_Real & sol);

#endif

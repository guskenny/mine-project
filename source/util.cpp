/* LP Relaxation: Open Pit Mining(PCPSP)
 * author: Davaa Baatar
 * date: 19.06.2016
*/
// utility functions
#include <util.h>
// verify solution
bool verify(const Daten & data, const Sol_Int & sol){
	const int r_max =data.getnResources();
	const int nBlocks = data.getNBlock();
	const int t_max = data.getNPeriod();
	const int d_max = data.getnDestination();
	const double rate = data.getDiscountRate();
	bool error = false;
// verify precedence constraints
	for(int a=0; a<nBlocks; a++){
		if(sol.x[a]<t_max){
			const Block & block_a = data.getBlock(a); //successor
			const std::vector<int> & preds = block_a.getPreds();
			const int n = preds.size();
			for(int p=0; p<n; p++){
				const int b = preds[p];    // predecessor
				if(sol.x[a] < sol.x[b] && sol.x[b] != t_max){
						std::cout<<(boost::format("ERROR: blocks %i and %i violated precedence constraint") %a %b )<<std::endl;
						std::cout<<(boost::format("\t pred %i is in period %i and succ. %i is in period %i") %b %sol.x[b] %a %sol.x[a])<<std::endl;
						error = true;
						return error;
				}
			}
			if(error) break;
		}
	}

// // verify distribution
// 	if(error) return error;
// 	for(int b=0; b<nBlocks; b++){
// 		if(sol.x[b] > -1 && sol.x[b] < data.getNPeriod()){
// 			double amount = 0.0;
// 			for(int d=0; d<d_max; d++) amount += sol.y[b][d];
// 			if(!(-eps_tol < amount - 1 && amount - 1 < eps_tol)){
// 					std::cout<<(boost::format("ERROR: distribution does not match!")) <<std::endl;
// 					std::cout<<(boost::format("\t block %i distributed %14.6f ") %b %amount) <<std::endl
// 									 <<"sol.x[" << b << "] = " << sol.x[b] << std::endl;
// 					error = 1;
// 					break;
// 			}
// 		}
// 	}
// verify not more than 1
	/// can not be violated for any integer solution

// verify resource constraints
	if(error) return error;
	for(int r=0; r<r_max; r++){
		long double total_use = 0;
		long double total_available = 0;
		for(int t=0; t<t_max; t++){
			long double amount = 0.0;

			for(int b=0; b<nBlocks; b++)
				if(sol.x[b]==t){
					const Block & block = data.getBlock(b);
					long double coef = block.getRCoef(0,r);
					if(fabs(coef) > 1e-20){
            amount += coef;
          }
			}
			// only for type L
			// get limits
			double R = data.getLimit(r, t);
			total_use += amount;
			total_available += R;
			// std::cout << "Period " << t << " resource (" << r << ") usage: " << amount << "/" << R << std::endl;
			// ckeck constraints
			if(amount > R + res_eps_tol){
				std::cout<<(boost::format("ERROR: Resource constraint is violated!")) <<std::endl;
				std::cout<<(boost::format(" resource r = %i t = %i   amount %14.6f    limit = %14.6f ") %r %t %amount %R)<<std::endl;
				error = true;
				return error;
			}
	}
	// std::cout << "\nTotal resource (" << r << ") usage: " << total_use << "/" << total_available << std::endl << std::endl;
	if(error) break;
}

// verify number of periods
	// if(error) return error;
	// int nPeriod = 0;
	// std::vector<bool> count(t_max,0);
	// for(int b=0; b<nBlocks; b++)
	// 	if(sol.x[b] > -eps_tol){
	// 		if(!count[sol.x[b]]){
	// 			count[sol.x[b]] = 1;
	// 			 nPeriod++;
	// 		}
	// 	}
	//
	//
	// if(nPeriod!=sol.nT){
	// 			std::cout<<(boost::format("WARNING: Number of periods do not match!")) <<std::endl;
	// 			std::cout<<(boost::format(" period in sol = %i actual number of periods  %i") %sol.nT %nPeriod)<<std::endl;
	// }

// verify objective value
	// if(error) return error;
	double factor = 1.0;
	double amount = 0.0;
	for(int t=0; t<t_max; t++){
		for(int b=0; b<nBlocks; b++){
			if(sol.x[b] == t){
				//ddfor(int d=0; d<d_max; d++)
				amount += data.getBlock(b).getProfit(0)*factor;
			}
		}
		factor /= (1+rate);
	}
	if(!(-obj_eps_tol < amount - sol.obj && amount - sol.obj <obj_eps_tol) ) {
				std::cout<<(boost::format("WARNING: Objective values do not match!" ))<<std::endl;
				std::cout<<(boost::format(" obj in sol = %14.6f  actual amount %14.6f") %sol.obj %amount) <<std::endl;
	}
return error;
}
bool verify(const Daten & data, const Sol_Real & sol){
	const int r_max =data.getnResources();
	const int nBlocks = data.getNBlock();
	const int t_max = data.getNPeriod();
	const int d_max = data.getnDestination();
	const double rate = data.getDiscountRate();
	bool error = 0;
// verify precedence constraints
	for(int a=0; a<nBlocks; a++){
		const Block & block_a = data.getBlock(a);  //succ.
		const std::vector<int> & preds = block_a.getPreds();
		const int n = preds.size();
		for(int p=0; p<n; p++){
			const int b = preds[p];
			const Block & block_b = data.getBlock(b); //pred
			double sum_a = 0.0;    //suc.
			double sum_b = 0.0;    //pred.
			for(int t=0; t<t_max; t++){
				sum_a += sol.x[a][t];
				sum_b += sol.x[b][t];
				if(sum_a > sum_b + eps_tol){
					std::cout<<(boost::format("ERROR: blocks %i and %i violated precedence constraint at period %i") %a %preds[p] %t)<<std::endl;
					std::cout<<(boost::format("\t sum of succ. and pred. resp. %14.6f %14.6f ") %sum_a %sum_b)<<std::endl;
					error = 1;
					break;
				}
			}
			if(error) break;
		}
		if(error) break;
	}

// verify distribution
	if(error) return error;
	for(int b=0; b<nBlocks; b++)
		for(int t=0; t<t_max; t++){
			double amount = 0.0;
			for(int d=0; d<d_max; d++) amount += sol.y[b][t][d];
			if(!(-eps_tol < amount - sol.x[b][t] && amount - sol.x[b][t] < eps_tol)){
				std::cout<<(boost::format("ERROR: distribution does not match!")) <<std::endl;
				std::cout<<(boost::format("\t period %i .x = %14.6f    but distributed %14.6f ") %t %sol.x[b][t] %amount) <<std::endl;
				error = 1;
				break;
			}
		if(error) break;
	}
// verify not more than 1
	if(error) return error;
	for(int b=0; b<nBlocks; b++){
		double amount = 0.0;
		for(int t=0; t<t_max; t++) amount += sol.x[b][t];
		if(amount > 1 + eps_tol){
				std::cout<<(boost::format("ERROR: More than 1 unit block utilized!")) <<std::endl;
				std::cout<<(boost::format("  block %i    amount %14.6f" ) %b %amount) <<std::endl;
				error = 1;
				break;
		}
	}

// verify resource constraints
	if(error) return error;
	for(int r=0; r<r_max; r++){
		for(int t=0; t<t_max; t++){
			double amount = 0.0;
			for(int b=0; b<nBlocks; b++){
				const Block & block = data.getBlock(b);
				for(int d=0; d<d_max; d++){
					amount += sol.y[b][t][d]*block.getRCoef(d, r);
				}
			}
			// only for type L
			// get limits
			double R = data.getLimit(r, t);
			// ckeck constraints
			if(amount > R + eps_tol){
				std::cout<<(boost::format("ERROR: Resource constraint is violated!")) <<std::endl;
				std::cout<<(boost::format(" resourc r = %i t = %i   amount %14.6f    limit = %14.6f ") %r %t %amount %R)<<std::endl;
				error = 1;
				break;
			}
	}
	if(error) break;
}

// verify number of periods
	if(error) return error;
	int nPeriod = 0;
	for(int t=0; t<t_max; t++)
		for(int b=0; b<nBlocks; b++){
			if(sol.x[b][t]>eps_tol){
				nPeriod++;
				break;
			}
	}
	if(nPeriod!=sol.nT){
				std::cout<<(boost::format("WARNING: Number of periods do not match!")) <<std::endl;
				std::cout<<(boost::format(" period in sol = %i actual number of periods  %i") %sol.nT %nPeriod)<<std::endl;
	}

// verify objective value
	if(error) return error;
	double factor = 1.0;
	double amount = 0.0;
	for(int t=0; t<t_max; t++){
		for(int b=0; b<nBlocks; b++)
			for(int d=0; d<d_max; d++){
				amount += data.getBlock(b).getProfit(d)*factor *sol.y[b][t][d];
		}
		factor /= (1+rate);
	}
	if(!(-eps_tol < amount - sol.obj && amount - sol.obj <eps_tol) ) {
				std::cout<<(boost::format("WARNING: Objective values do not match!" ))<<std::endl;
				std::cout<<(boost::format(" obj in sol = %14.6f  actual amount %14.6f") %sol.obj %amount) <<std::endl;
	}
return error;
}
int convert_triplet_to_index(const Daten & data, int b, int t, int d){
	/// index = b*t_max*d_max + t*t_max + d
	const int t_max = data.getNPeriod();
	const int d_max = data.getnDestination();
	int index = b*t_max*d_max + t*d_max + d;
	return index;
}
void convert_index_to_triplet(const Daten & data, int index, int & b, int & t, int & d){
	const int t_max = data.getNPeriod();
	const int d_max = data.getnDestination();
	b = (int) index/(t_max*d_max);
	index -= b*t_max*d_max;
	t = (int) index/d_max;
	index -= t*d_max;
	d = index;
}
void convert_Z_to_XY(Daten & data, std::vector<double> & sol_z, std::vector<std::vector<double> > & X, \
					std::vector< std::vector< std::vector<double> > > & Y, int & nPeriods){
	const int nBlocks = data.getNBlock();
	const int t_max = data.getNPeriod();
	const int d_max = data.getnDestination();
	int n = nBlocks*t_max*d_max;

	// get Y values
	for(int i=0; i<n; i++){
		int b,t,d;
		//int zb, zt, zd;
		convert_index_to_triplet(data, i, b, t, d);
		if(d>0){
			int index_1 = convert_triplet_to_index(data, b, t, d);
			int index_2 = convert_triplet_to_index(data, b, t, d-1);
			Y[b][t][d] = sol_z[index_1] - sol_z[index_2];
		}else if (t>0){
				 int index_1 = convert_triplet_to_index(data, b, t, 0);
				 int index_2 = convert_triplet_to_index(data, b, t-1, d_max-1);
				 Y[b][t][d] = sol_z[index_1] - sol_z[index_2];
			 }else {
				 int index = convert_triplet_to_index(data, b, 0, 0);
				 Y[b][t][d] = sol_z[index];
			 }
	}

	// get X values and number of periods
	nPeriods = 0;
	for(int b=0; b<nBlocks; b++)
		for(int t =0; t<t_max; t++){
			X[b][t] = 0.0;
			for(int d=0; d<d_max; d++)
				X[b][t] += Y[b][t][d];
//			if(t>=nPeriods && X[b][t] >= eps_tol) nPeriods++;
	}
	// push forward
	nPeriods = 0;
	for(int t = 0; t < t_max; t++){
		bool active = 0;
		for(int b = 0; b < nBlocks; b++)
			if(X[b][t] > eps_tol) {
				if(nPeriods != t){
					X[b][nPeriods] = X[b][t];
					X[b][t] = 0.0;
					for(int d=0; d<d_max; d++){
						Y[b][nPeriods][d] = Y[b][t][d];
						Y[b][t][d] = 0.0;
					}
				}
				if(!active) active = 1;
			}
		if(active) nPeriods++;
	}
}
double evaluateObjValue(const Daten & data, Sol_Real & sol){
	double factor = 1.0;
	double amount = 0.0;
	const int nBlocks = data.getNBlock();
	const int t_max = data.getNPeriod();
	const int d_max = data.getnDestination();
	const double rate = data.getDiscountRate();

	for(int t=0; t<t_max; t++){
		for(int b=0; b<nBlocks; b++)
			for(int d=0; d<d_max; d++){
				amount += data.getBlock(b).getProfit(d)*factor *sol.y[b][t][d];
		}
		factor /= (1+rate);
	}
	sol.obj = amount;
	return amount;
}
void info_data(Daten data, char * filename){

	const char pT = data.getProbType();
	const int nB = data.getNBlock();
	const int t_max = data.getNPeriod();
	const int d_max = data.getnDestination();
	const int r_max = data.getnResources();
	const double rate = data.getDiscountRate();

	// info on screen
	std::string datafile(filename);
	std::string strng;
	std::string::size_type foundPos, foundEnd;
	foundPos = datafile.find_last_of("/");
	if(foundPos==std::string::npos) strng="";
	else strng = datafile.substr(foundPos+1);
	printf("----------------------------------------------------------------------------------------------------------\n");
	printf("Data: %s\n", strng.c_str());
	printf("Type: ");
	if(pT=='p') printf("PCPSP\n");
	else if (pT=='c') printf("CPIP\n");
		 else printf("UPIP\n");
	printf("NBlocks: %i\n", nB);
	if(pT!='u'){
		printf("NPeriods: %i\n", t_max);
		printf("NResource_Side_Constraints: %i\n",r_max);
		printf("Discount_Rate: %.6f \n",rate);
		if(pT!='c'){
			printf("NDestination: %i\n",d_max);
		}
	}
}
void usage(char ** argv ){
	 std::cerr << argv[0] << " [options] Path\n"
	 << "Solve PCPSP (mine planning optimisation) using a cummulative variable formulation\n"
       << "Path      directory containing data files (uncompressed)\n"
       << "-C        CPLEX solver (not implemented yet)\n"
       << "-G	     Gurobi solver (not implemented yet)\n"
       << "-V	     Volume method (not implemented yet)\n"
       << "-L	     Lagrangian optimisation (not implemented yet)\n"
       << "-T	     Lagrangian optimisation \n"
       << "-n	     newtork flow max closure solver\n"
       << "-b	     Boykov-Kolmogorov max flow algorithm for max closure\n"
       << "-e	     Edmonds-Karp max flow algorithm\n"
       << "-p	     Push-relable max flow algorithm\n"
       << "-w fn     Write cumulative model to file 'fn'\n"
       << std::endl;
}
double evaluate_const(Daten &data, std::vector< std::vector< double > >  &lambda){
	std::vector< std::vector< std::vector<int> > > * limits = data.getLimit();
	double obj_const = 0.0;
	const int r_max =data.getnResources();
	const int t_max = data.getNPeriod();
	for(int r=0; r<r_max; r++)
		for(int t=0; t<t_max; t++)
			obj_const += lambda[t][r]*(*limits)[r][t][0];
	return obj_const;
}
void evaluate_supply(Daten & data, std::vector<std::vector<double> > & lambda, double * supply){
	const int r_max =data.getnResources();
	const int nBlocks = data.getNBlock();
	const int t_max = data.getNPeriod();
	const int d_max = data.getnDestination();
	const double rate = data.getDiscountRate();
	std::vector<Block> * blocks = data.getBlock();
	double total_demand = 0;
	for(int b=0; b<nBlocks; b++){
		std::vector<double> profit((*blocks)[b].getProfit()->begin(), (*blocks)[b].getProfit()->end()); // copy the profit
		std::vector< std::vector<double> >  * q = (*blocks)[b].getRCoef();
		for(int t=0; t<t_max; t++){
			for(int d=0; d<d_max; d++){
				int index = convert_triplet_to_index(data, b,t,d)+1;
				supply[index] = profit[d];

				if(d==d_max-1 && t<t_max-1)
					supply[index] -= profit[0]/(1+rate);
				else if(d<d_max-1)
						supply[index] -= profit[d+1];

				if(!(d==d_max-1 && t==t_max-1)){ // negative part for L
					for(int r=0; r<r_max; r++)
						supply[index] -= lambda[t][r]*((*q)[d][r]);
				}

				if(d==d_max-1 && t<t_max-1) { // positive part for L
					for(int r=0; r<r_max; r++){
						supply[index] += lambda[t+1][r]*((*q)[0][r]);}
				}else if(d!=d_max-1 && !(d==0 && t==t_max-1)){
						for(int r=0; r<r_max; r++)
							supply[index] += lambda[t][r]*((*q)[d+1][r]);
				}
				total_demand += supply[index];
				supply[index] *= -1;
			}
			for(int d=0; d<d_max; d++) profit[d] /= (1+rate); // prepare next profit rate
		}
	}
	supply[0] = total_demand;    // supply of the dummy node
}

void save_sol(char * filename, char opt, bool status, const double obj_val,
	const double ub, const std::vector<std::vector<double> > & X,
	const std::vector< std::vector< std::vector<double> > > & Y, const int nPeriods,
	const int nIter, const double total_time, const double total_MC_time,
	const double cpu_network, const double total_RLP_time, std::string attr, std::string text){
	// output save solution
	const int nB = (int) X.size();
	const int t_max = (int) X[0].size();

	const int d_max = (int) Y[0][0].size();

	 std::string datafile(filename);
	 std::string strng;
	 std::string::size_type foundPos;
	 foundPos = datafile.find_last_of("/");
	 if(foundPos==std::string::npos) strng="";
	 else strng = datafile.substr(foundPos+1);
	 strng = strng + "_LPR_BZ.sol";


	 std::ofstream fp_out;
	 fp_out.open(strng.c_str(), std::ios::app);
	 if(fp_out.is_open()){
		 fp_out<< text <<std::endl;
		 fp_out<<"Max Closure Solver: ";
			switch (opt) {
			  case 'n': case 'N':
					fp_out<<"Network Simplex Algorithm "<< std::endl;
					break;
			  case 'b': case 'B':
					fp_out<<"Boykov-Kolmogorov max flow algorithm "<< std::endl;
					break;
			  case 'e': case 'E':
					fp_out<<"Edmonds-Karp max flow algorithm "<< std::endl;
					break;
			  case 'p': case 'P':
					fp_out<<"Push-relable max flow algorithm "<< std::endl;
					break;
			  default:
				  std::cerr << "Unknown option '" << (char) opt << "'\n";
			}
		 fp_out<<(boost::format("\n%-36s %s") %"Status" %( (status) ? "Optimal " : "Feasible/NotSolved" ))<<std::endl;
		 fp_out<<(boost::format("%-36s %i") %"Number of iterations" %nIter)<<std::endl;
		 fp_out<<(boost::format("%-36s %-8.2f sec.") %"Total CPU time" %total_time)<<std::endl;
		 fp_out<<(boost::format("   %-33s %8.2f sec.") %"Network is constructed in" %cpu_network)<<std::endl;
		 fp_out<<(boost::format("   %-33s %8.2f sec.") %"Max Closure: Total CPU time" %total_MC_time)<<std::endl;
		 fp_out<<(boost::format("   %-33s %8.2f sec.") %"Restricted LP: Total CPU time" %total_RLP_time)<<std::endl;

		 fp_out<<(boost::format("%-36s %-28.6f") %"Objective value" %obj_val)<<std::endl;
		 fp_out<<(boost::format("   %-33s %-9.6f %s") %"Gap" %((ub - obj_val)/obj_val) %"%")<<std::endl;
		 fp_out<<(boost::format("\n%-36s %i") %"Number of periods" %nPeriods)<<std::endl;
		 /* solution file format blocks -> time period -> destinations */
		 if(status){
			 fp_out<<(boost::format("\n%=6s  %=8s  %=2s  %12s")%"Block" %"X" %"t" %"destinations")<<std::endl;
			 for(int b=0; b<nB; b++){
				fp_out<<(boost::format("%6i") %b);
				bool YesNo = 0;
				for(int t = 0; t<t_max; t++)
					if(X[b][t]> eps_tol){
						if(YesNo) fp_out <<(boost::format("      "));
						fp_out<<(boost::format("  %8.6f  %2i")%X[b][t] %t);
						for(int d = 0; d<d_max; d++){
							fp_out<<(boost::format("  %8.6f")%Y[b][t][d]);
						}
						fp_out<<std::endl;
					YesNo = 1;
				// break;
				}
				if(!YesNo) fp_out<<(boost::format("  %8i")%0)<< std::endl;
			//	else fp_out<<std::endl;
			 }
		}
		 fp_out.close();
	 }
}
void save_sol_int(char * filename, bool status, const Sol_Int & sol, const double obj_val, const double total_time, std::string attr, std::string text){
		// output save solution
	const int nB = (int) sol.x.size();
	const int d_max = (int) sol.y[0].size();

	 std::string datafile(filename);
	 std::string strng;
	 std::string::size_type foundPos;
	 foundPos = datafile.find_last_of("/");
	 if(foundPos==std::string::npos) strng="";
	 else strng = datafile.substr(foundPos+1);
	 strng = strng + attr;

	 std::ofstream fp_out;
	 fp_out.open(strng.c_str(), std::ios::app);
	 if(fp_out.is_open()){
		 fp_out<< text <<std::endl;
		 fp_out<<(boost::format("\n%-36s %s") %"Status" %( (status) ? "Optimal " : "Feasible/NotSolved" ))<<std::endl;
		 fp_out<<(boost::format("%-36s %-8.2f sec.") %"Total CPU time" %total_time)<<std::endl;

		 fp_out<<(boost::format("%-36s %-28.6f") %"Objective value" %sol.obj)<<std::endl;
		 fp_out<<(boost::format("   %-33s %-9.6f") %"Gap" %((obj_val - sol.obj)/sol.obj) )<<std::endl;

		 fp_out<<(boost::format("\n%-36s %i") %"Number of periods" %sol.nT)<<std::endl;
		 /* solution file format blocks -> time period -> destinations */
		 if(status){
			 fp_out<<(boost::format("\n%=6s  %=8s  %=2s  %12s")%"Block" %"X" %"t" %"destinations")<<std::endl;
			 for(int b=0; b<nB; b++){
				fp_out<<(boost::format("%6i") %b);
				if(sol.x[b]> -eps_tol){
						fp_out<<(boost::format("  %8i  %2i")%1 %sol.x[b]);
						for(int d = 0; d<d_max; d++){
							fp_out<<(boost::format("  %8.6f")%sol.y[b][d]);
						}
						fp_out<<std::endl;
				}else{
					fp_out<<(boost::format("  %8i")%0)<< std::endl;
				}
			}
		}
		 fp_out.close();
	 }
}

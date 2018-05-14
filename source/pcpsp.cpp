/***************************************************************
 * 
 * MIP model: Precedence Constrained Product Scheduling Problem
 * 
 * ***********************************************************/
 #include <ilcplex/ilocplex.h>
 #include <daten.h>
 #include <cmath>
 //#include <math.h> 
 ILOSTLBEGIN
 typedef IloArray< IloBoolVarArray > IloBoolVarArray2;
 typedef IloArray< IloNumVarArray > IloNumVarArray2;
 typedef IloArray< IloNumVarArray2 > IloNumVarArray3;
int main(int argc, const char** argv){
 IloEnv env;
 IloTimer timer(env);
 if(argc < 1){
	std::cout<<"Usage : ./exe [OPTION] PATH "<<std::endl;
	std::cout<<"       Options"<<std::endl;
	std::cout<<"          -p    PCPSP"<<std::endl;
	std::cout<<"Sample "<<std::endl;
	std::cout<<"      ./mining -p Data/zuck_small"<<std::endl;
	exit(-1);
 }
 if(argv[1][1]!='p'){
	 cout<<"Unknown Problem Type!"<<endl;
	 exit(-1);
 }
 try{
	Daten data(argv[2], argv[1][1]);
//	Daten data(argv[2]);
	const char pT = data.getProbType();
	const int nB = data.getNBlock();
	const int t_max = data.getNPeriod();
	const int d_max = data.getnDestination();
	const int r_max = data.getnResources();
	const double rate = data.getDiscountRate();
	std::vector<Block> * blocks=data.getBlock();
	cout<<setprecision(12);
	 // info on screen
	 {
		std::string datafile(argv[2]);
		std::string strng;
		std::string::size_type foundPos, foundEnd;
		foundPos = datafile.find_last_of("/");
		if(foundPos==std::string::npos) strng="";
		else strng = datafile.substr(foundPos+1);
		cout<<"Data: " << strng << endl;
		cout<<"Type: ";
		if(pT=='p') cout<< "PCPSP"<<endl;
		else if (pT=='c') cout<< "CPIP"<<endl;
			 else cout<< "UPIP"<<endl;
		cout<<"NBlocks: "<< nB<<endl;
		if(pT!='u'){
			cout<<"NPeriods: "<< t_max<<endl;
			cout<<"NResource_Side_Constraints: "<< r_max<<endl;
			cout<<"Discount_Rate: "<< rate<<endl;	 
			if(pT!='c'){
				cout<<"NDestination: "<< d_max<<endl;
			}
		}
	 }
	 // decision variables
	 IloBoolVarArray2 Xvar(env); // block -> period => usage
	 IloNumVarArray3 Yvar(env);  // block -> period ->destination => value(%)
	 for(int b=0; b<nB; b++){
		 Xvar.add(IloBoolVarArray(env, t_max));
	 }
	 for(int b=0; b<nB; b++){
		 Yvar.add(IloNumVarArray2(env));
		 for(int t=0; t<t_max; t++)
			Yvar[b].add(IloNumVarArray(env, d_max, 0, 1, ILOFLOAT));
	 }

	 // type of constraints
	 /*Precedence (7) */
	 IloRangeArray Precedence(env);
	 for(int a=0; a<nB; a++){
		 std::vector<int> * pred = (*blocks)[a].getPreds();
		 int n = (*blocks)[a].getNumPred();
		 for(int p=0; p<n; p++){
			 int b = (*pred)[p];
			 for(int t=0; t<t_max; t++){
				 IloExpr expr(env);
				 for(int j=0; j<=t; j++)
					expr += Xvar[b][j] - Xvar[a][j];
				 Precedence.add(IloRange(env, 0, expr, 1));
			 }
		}
	 }
	
	 /* SumDest (8) */
	 IloRangeArray SumDest(env);
	 for(int b=0; b<nB; b++){
		for(int t=0; t<t_max; t++){
			IloExpr expr(env);
			expr -= Xvar[b][t];
			for(int d=0; d<d_max; d++)
				expr += Yvar[b][t][d];
			 SumDest.add(IloRange(env, 0, expr, 0));
		 }
	 }
	 
	 /* CliqueBlock (9) */
	 IloRangeArray CliqueBlock(env);
	 for(int b=0; b<nB; b++){
		 IloExpr expr(env);
		 for(int t=0; t<t_max; t++)
			expr += Xvar[b][t];
		 CliqueBlock.add(IloRange(env, 0, expr, 1));	
	 }
	 /* Resource constraints (10) */
	 IloRangeArray Resource(env);
	 for(int r=0; r<r_max; r++)
		for(int t=0; t<t_max; t++){
			char cType = data.getResConstrType(r, t);
			IloExpr expr(env);
			for(int b=0; b<nB; b++){
				for(int d=0; d<d_max; d++){
					double coef = (*blocks)[b].getRCoef(d,r);
					expr += coef*Yvar[b][t][d];
				}
			}
			if(cType == 'L'){ 
				Resource.add(IloRange(env, -IloInfinity, expr, data.getLimit(r, t)));
			}else if(cType == 'R'){
					Resource.add(IloRange(env, data.getLimit(r, t), expr, IloInfinity));
				 }else if(cType == 'I'){
					 Resource.add(IloRange(env, data.getLimit(r, t), expr, data.getLimit(r, t, 1)));
				 }
	 }
	 
	 // objective function
	 IloExpr exprObj(env);
	 for(int b=0; b<nB; b++){
		for(int d=0; d<d_max; d++){
			double coef = (*blocks)[b].getProfit(d);
			for(int t=0; t<t_max; t++){
				exprObj += coef*Yvar[b][t][d];
				coef /= (1+rate); 
			}
		}
	 }
	 IloObjective Obj=IloMaximize(env, exprObj);
	 // cplex model 
	 IloModel pcpsp(env);
	 pcpsp.add(Precedence);
	 pcpsp.add(SumDest);
	 pcpsp.add(CliqueBlock);
	 pcpsp.add(Resource);
	 pcpsp.add(Obj);
	 IloCplex pcpsp_cpx(pcpsp);
	 pcpsp_cpx.setParam(IloCplex::ClockType,1);
	 
	 // solve & retrieve related info
	 timer.start();
	 pcpsp_cpx.solve();
	 timer.stop();
	 // output save solution
	 
	 std::string datafile(argv[2]);
	 std::string strng;
	 std::string::size_type foundPos;
	 foundPos = datafile.find_last_of("/");
	 if(foundPos==std::string::npos) strng="";
	 else strng = datafile.substr(foundPos+1);
	 strng = strng + "_pcpsp.sol";
	 
	 ofstream fp_out;
	 fp_out.open(strng.c_str(), ios::app);
	 if(fp_out.is_open()){
		 fp_out<< "Status          "<<pcpsp_cpx.getStatus()<< endl;
		 fp_out<< "CPU time        "<<timer.getTime()<<" sec."<< endl;
		 fp_out<< "Objective Value "<<pcpsp_cpx.getObjValue() << endl;
		 fp_out<< "Block X t destinations"<<endl;
		 /* solution file format blocks -> time period -> destinations */
		 for(int b=0; b<nB; b++){
			fp_out<<b;
			IloNumArray sol(env);
			pcpsp_cpx.getValues(sol, Xvar[b]);
			bool YesNo = 0;
			for(int t =0; t<t_max; t++)
				if(sol[t]>0.9){
					fp_out<<" "<< 1 << " "<< t ;
					IloNumArray distr(env);
					pcpsp_cpx.getValues(distr, Yvar[b][t]);	
					for(int d =0; d<d_max; d++){
						fp_out<<" "<< distr[d];
					}
				YesNo = 1;
				break;
			}
			if(!YesNo) fp_out<<" "<< 0 << endl;
			else fp_out<<endl;
		 }
		 fp_out.close();	 
	 }
	// clean
 }
 catch (IloException & ex) { cerr << "Error: " << ex << endl; }
 catch (...) { cerr << "Unknown Error Occured" << endl;}
 env.end();
 
 return 0;
}

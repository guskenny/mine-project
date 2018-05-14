/***************************************************************
 *
 * MIP model: Precedence Constrained Product Scheduling Problem
 *
 * ***********************************************************/
#include <boost/format.hpp>
#include "QOL/QolMIP.h"
#include "QOL/QolColFormulation.h"
#include "QOL/CplexFormulation.h"
#include "QOL/GurobiFormulation.h"
//#include "QOL/CPBeamACO.h" // needs compilation with boost::thread
#include "QOL/CpuTimer.h"
#include <daten.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <SinglePSolver.h> // remove to return to original
//#include <math.h>
int main(int argc, const char** argv){
//  qol::MIPSolver *mipPtr=0; // uncomment for original


  if(argc <= 1){
    std::cout<<"Usage : ./QOL [OPTION] PATH "<<std::endl;
    std::cout<<"       Options"<<std::endl;
    //std::cout<<"          -p    PCPSP"<<std::endl;
    std::cout<<"Sample "<<std::endl;
    std::cout<<"      ./mining -g Data/zuck_small"<<std::endl;
    std::cout<<"Options: g/c solve with gurobi/cplex\n";
    std::cout<<"-r solve relaxed LP only\n";
    exit(-1);
  }
  try{
    SinglePSolver solver = SinglePSolver(argc,argv);
    solver.solve();

/* COMMENTED ORIGINAL CODE OUT TO TEST SINGLE PERIOD MODEL
 * REMOVE TWO LINES ABOVE AND UNCOMMENT EVERYTHING BELOW
 * AND LINE JUST BELOW MAIN DEFINITION TO RETURN TO ORIGINAL */

    //	Daten data(argv[2], argv[1][1]);
    Daten data(argv[argc-1],'p');
    const char pT = data.getProbType();
    const int nB = data.getNBlock();
    const int t_max = data.getNPeriod();
    const int d_max = data.getnDestination();
    const int r_max = data.getnResources();
    const double rate = data.getDiscountRate();
    std::vector<Block> * blocks=data.getBlock();
    std::cout<<std::setprecision(12);
    // info on screen
    {
      std::string datafile(argv[argc-1]);
      std::string strng;
      std::string::size_type foundPos, foundEnd;
      foundPos = datafile.find_last_of("/");
      if(foundPos==std::string::npos) strng="";
      else strng = datafile.substr(foundPos+1);
      std::cout<<"Data: " << strng << std::endl;
      std::cout<<"Type: ";
      if(pT=='p') std::cout<< "PCPSP"<<std::endl;
      else if (pT=='c') std::cout<< "CPIP"<<std::endl;
      else std::cout<< "UPIP"<<std::endl;
      std::cout<<"NBlocks: "<< nB<<std::endl;
      if(pT!='u'){
	std::cout<<"NPeriods: "<< t_max<<std::endl;
	std::cout<<"NResource_Side_Constraints: "<< r_max<<std::endl;
	std::cout<<"Discount_Rate: "<< rate<<std::endl;
	if(pT!='c'){
	  std::cout<<"NDestination: "<< d_max<<std::endl;
	}
      }
    }
    bool solveRelaxed= false;
    std::string lpFile="",solnFile=data.getName()+".sol";
    qol::Parameters param;
    param.setParamVal(qol::VERBOSITY,1);
    for(int c=1;c<argc-1;++c){
      std::string flag=argv[c];
      if(flag == "-g")
	mipPtr=new qol::GurobiFormulation();
      else if(flag=="-c")
	mipPtr=new qol::CplexFormulation();
      // else if(flag=="-b")
      //	mipPtr=new qol::CPBeamACO();
      else if(flag=="-r")
	solveRelaxed = true;
      else if(flag=="-w" )
	lpFile=std::string(argv[c+1]);
      else if(flag=="-o")
	solnFile=std::string(argv[c+1]);
      else if(flag =="-t")
	param.setParamVal(qol::TIMELIMIT,atof(argv[++c]));
      else
	std::cerr<<"WARNING: unknown flag " << flag << " ignored\n";
    }
    if(! mipPtr) {
      std::cout << "No solver specified - use flags -g|-c|/-b\n";
      exit(2);
    }
    qol::MIPSolver &mip=*mipPtr;
    mip.setParameters(param);
    // decision variables
    mip.setVarSize(nB*t_max*(1+d_max)); // mainly for efficiency
    std::vector<std::vector<qol::Variable> > x(nB);
    std::vector<std::vector<std::vector<qol::Variable> > > y(nB);


    for(int b=0; b<nB; b++){
      x[b].resize(t_max);
      for(int t=0; t<t_max; t++){
	x[b][t]=mip.addVar(0,1,0.0,qol::Variable::BINARY, // LB,UB,cost
			   boost::str(boost::format("x%03d_%03d")%b%t));
      }
    }
    for(int b=0; b<nB; b++){
      y[b].resize(d_max);
      for(int d=0;d<d_max;++d){
	y[b][d].resize(t_max);
	double profit=(*blocks)[b].getProfit(d);
	for(int t=0; t<t_max; t++){
	  y[b][d][t]=mip.addVar(0,1,-profit,qol::Variable::CONTINUOUS,
				boost::str(boost::format("y%d_%d_%d")%b%d%t));
	  // could write objective as sum expression but this is faster and easy
	  profit /= (1+rate);
	}
      }
    }


    // type of constraints
    /*Precedence (7) (but using different definition of x */
    for(int a=0; a<nB; a++){
      std::vector<int> * pred = (*blocks)[a].getPreds();
      int n = (*blocks)[a].getNumPred();
      for(int p=0; p<n; p++){
	int b = (*pred)[p];
	for(int t=0; t<t_max; t++)
	  mip.addConstraint(x[b][t] >= x[a][t]).setName("pre%d_%d_%d",b,a,t);
      }
    }

    /* SumDest (8) */
    for(int b=0; b<nB; b++){
      for(int t=0; t<t_max; t++){
	qol::Expression sumY;
	for(int d=0; d<d_max; d++)
	  sumY += y[b][d][t];
	qol::Expression rhs=x[b][t];
	if(t > 0) rhs -= x[b][t-1];
	mip.addConstraint( sumY == rhs ).setName("dest%d_%d",b,t);
      }
    }

    /* Block remains done (9) */
    for(int b=0; b<nB; b++){
      for(int t=1; t<t_max; t++)
	mip.addConstraint( x[b][t-1] <= x[b][t]).setName("xt%d_%d",b,t);
    }
    /* Resource constraints (10) */
    for(int r=0; r<r_max; r++)
      for(int t=0; t<t_max; t++){
	char cType = data.getResConstrType(r, t);
	qol::Expression expr;
	for(int b=0; b<nB; b++){
	  for(int d=0; d<d_max; d++){
	    double coef = (*blocks)[b].getRCoef(d,r);
	    if(fabs(coef) > 1e-5) expr += coef*y[b][d][t];
	  }
	}
	if(cType == 'L'){
	  mip.addConstraint(expr <= data.getLimit(r, t)).setName("R%d_%d",r,t);
	}else if(cType == 'G'){
	  mip.addConstraint(expr >= data.getLimit(r, t)).setName("R%d_%d",r,t);
	}else{
	  std::cerr << "ERROR: resource constraint type " << cType << " not implemented - IGNORED\n";
	}
      }

    if(lpFile != ""){
      mip.writeLP(lpFile.c_str());
      std::cout << "Wrote " << lpFile << std::endl;
    }
    qol::CpuTimer timer;
    qol::Status status = solveRelaxed ? mip.solveRelaxed() : mip.solveExact();
    std::cout << boost::format("Completed in %.2f sec CPU / %.2f sec wall. Objective = %f\n"
			       ) % timer.elapsedSeconds() % timer.elapsedWallTime() % (-mip.getObjective());

    std::ofstream fp_out;
    fp_out.open(solnFile.c_str()); //, std::ios::app);
    if(fp_out.is_open()){
      std::cout << "Writing solution to " << solnFile << std::endl;
      fp_out<< "# Status          "<<status<< std::endl;
      fp_out<< "# CPU time        "<<timer.elapsedSeconds()<<" sec.  Wall:"<< timer.elapsedWallTime() <<std::endl;
      fp_out<< "# Objective Value "<<-mip.getObjective() << std::endl;
      fp_out<< "# Block destination time yval"<<std::endl;
      /* solution file format blocks -> destinations -> time non-zero y value*/
      for(int b=0; b<nB; b++){
	for(int t=0;t<t_max;++t)
	  for(int d=0;d<d_max;++d){
	    if(mip.getPrimal(y[b][d][t]) > 1e-5)
	      fp_out << b << " " << d << " " << t << " " << mip.getPrimal(y[b][d][t]) << std::endl;
	  }
      }
      fp_out.close();
    }
    // clean
  } // end try statement
  catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
  catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}
  return 0;
}

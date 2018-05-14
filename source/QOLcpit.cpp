/***************************************************************
 *
 * MIP model: CPIT
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

  } // end try statement
  catch (qol::Exception & ex) { std::cerr << "Error: " << ex.what() << std::endl; }
  catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}
  return 0;
}

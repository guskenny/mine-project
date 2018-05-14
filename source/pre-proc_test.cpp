#include "network.h"
#include "../include/boostMaxFlow.h"
#include "../include/MaxClosure_Base.h"
#include "../include/graph.h"
#include <map>
#include <ctime>
#include "../include/MaxClosure_PP.h"
#include <iostream>
#include <cstdlib>

int
main()
{
  // instantiate network 
  Network myNet = Network();

  // read DIMACS file to network
  std::vector<double> profits = myNet.readDimacsNoST();

  Graph myGraph = myNet.makeGraph(); 

  std::map<int,int> fixed;

  //fixed[37] = 1;
  //fixed[31] = 1;
  //fixed[44] = 0;

  //std::cout << "element 37 is " << fixed[37] << std::endl;

  MaxClosureFactory mcfactory;

  mcfactory.setOption('b');
  
  MaxClosure_PP pre_proc = MaxClosure_PP(myGraph,fixed,mcfactory);

  //  std::cout << endl << flowOO << endl;

  //maxFlow.printFlowOutput(flowOO);

  std::vector<int> myVec;

  myVec.push_back(3);
  myVec.push_back(1);
  myVec.push_back(7);
  myVec.push_back(4);

//  std::cout << "\nunsorted: ";
//  for (int i = 0; i < myVec.size(); i++)
//    std::cout << myVec[i] << " ";

//  std::cout << "\nsorted: ";
//  std::sort(myVec.rbegin(), myVec.rend());

//  for (int i = 0; i < myVec.size(); i++)
//    std::cout << myVec[i] << " ";

//  std::cout << std::endl;

  std::vector<double> profits2;
  float MIN = -3.0;
  float MAX = 6.0;
  srand(time(NULL));
 
  for (int i = 0; i < 100; i++){
    profits2.push_back(((MAX-MIN)*((float)rand()/RAND_MAX))+MIN);
  }

  pre_proc.setProfit(profits2);

  pre_proc.solve();

  std::vector<int> sol(100,0);

  pre_proc.getClosure(sol);

//  MaxClosure_Base reducedProblem = pre_proc.getReducedProblem();

//  int numReducedNodes = reducedProblem->getGraph()->getNumNodes();
  

  return 0;
}

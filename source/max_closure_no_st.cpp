/*************************************************************************** 
                          Main function to test boost max flow algorithm
                          ----------------------------------------------- 
    last modified   : 21/6/2016 
    copyright       : (C) 2016 by Angus Kenny
    libraries		    : . 
    description		  : tests boost max flow algorithm
 ***************************************************************************/ 
 
/*************************************************************************** 
 *                                                                         * 
 *   This program is free software; you can redistribute it and/or modify  * 
 *   it under the terms of the GNU General Public License as published by  * 
 *   the Free Software Foundation; either version 2 of the License, or     * 
 *   (at your option) any later version.                                   * 
 *                                                                         * 
 ***************************************************************************/ 

#include "network.h"
#include "../include/boostMaxFlow.h"
#include "../include/graph.h"
#include "../include/MaxClosure_BoostMaxFlow_EK.h"
#include "../include/MaxClosure_BoostMaxFlow_PR.h"
#include "../include/MaxClosure_BoostMaxFlow_BK.h"
#include "../include/MaxClosure_Base.h"

#include <ctime>

int
main()
{
  std::clock_t start,finish,total_start,total_finish,init_start,init_finish,build_start,build_finish;

  init_start = std::clock();

  // instantiate network 
  Network myNet = Network();

  // read DIMACS file to network
  std::vector<double> profits = myNet.readDimacsNoST();

  Graph myGraph = myNet.makeGraph(); 

  // instantiate BoostMaxFlow object
  BoostMaxFlow maxFlow = BoostMaxFlow();  

  // build boost network with generic network object
  maxFlow.buildNetworkNoST(&myNet, &profits);

  init_finish = std::clock();

  build_start = std::clock();

//  MaxClosure_BoostMaxFlow_EK myEK = MaxClosure_BoostMaxFlow_EK(myGraph);

//  myEK.setProfit(profits);
 
  MaxClosure_BoostMaxFlow_BK myBK = MaxClosure_BoostMaxFlow_BK(myGraph);

  myBK.setProfit(profits);
  
  build_finish = std::clock(); 

//  MaxClosure_BoostMaxFlow_PR myPR = MaxClosure_BoostMaxFlow_PR(myGraph);

//  myPR.setProfit(profits);
  
//  vector<int> closure = vector<int>(myGraph.getNumNodes(),0);

//  std::cout << "\n------------------------------\n"
//            << "Network width: " << myNet.getWidth()
//            << "\nNetwork depth: " << myNet.getDepth()
//            << "\nNumber of nodes: " << myNet.getNumNodes()
//            << "\nNumber of arcs: " << myNet.getNumArcs()
//            << "\n";

//  std::cout << "\nTime to read file: " << (init_finish-init_start) / (double)(CLOCKS_PER_SEC / 1000) << " ms\n";
//  std::cout << "Time to build graph: " << (build_finish-build_start) / (double)(CLOCKS_PER_SEC / 1000) << " ms\n\n";
  long flowOO;

  total_start = std::clock();

  start = std::clock();

  flowOO = myBK.solve();  

  finish = std::clock();

//  std::cout << "Boykov-Kolmogorov time: " << (finish-start) / (double)(CLOCKS_PER_SEC / 1000) << " ms\n";
 
//  start = std::clock();

//  flowOO = myPR.solve();  

//  finish = std::clock();

//  std::cout << "Push-relabel time: " << (finish-start) / (double)(CLOCKS_PER_SEC / 1000) << " ms\n";

//  start = std::clock();
//
//  flowOO = myEK.solve();  
//
//  finish = std::clock();
//
//  std::cout << "Edmonds-Karp time: " << (finish-start) / (double)(CLOCKS_PER_SEC / 1000) << " ms\n";

  start = std::clock();
  vector<int> closure = vector<int>(myNet.getNumNodes(),0);
  
  myBK.getClosure(closure);

  finish = std::clock();

  total_finish = std::clock();
  
//  std::cout << "\nTime to find closure: " << (finish-start) / (double)(CLOCKS_PER_SEC / 1000) << " ms\n";
//  std::cout << "\nTotal run time: " << (total_finish-total_start) / (double)(CLOCKS_PER_SEC / 1000) << " ms\n";

//  std::cout << "\n------------------------------\n";


//  // run push relabel algorithm
//  long flowOO = maxFlow.pushRelabel();
//

  // Call the recursive visiting function to print DFS traversal
  //maxFlow.getClosure(closure);

  maxFlow.printClosure(closure);

//  BoostMaxFlow::printClosureGraph(closure, 10,10);

//  std::cout << endl << flowOO << endl;

  //maxFlow.printFlowOutput(flowOO);

  return 0;
}

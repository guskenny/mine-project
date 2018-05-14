#include "UpitSolver.h"
#include "MaxClosureFactory.h"

UpitSolver::UpitSolver(const Daten &prob){
  nb = prob.getNBlock();

  //std::cout << "nb: " << nb << std::endl;

  // define the actual graph
  Graph upitGraph;
  upitGraph.resizeNodes(nb);

  // define profit vector
  std::vector<double> upitProfit;
  upitProfit.resize(nb);

  // iterate over all blocks
  for(int b=0; b < nb; ++b){
	  const Block &block = prob.getBlock(b);
    // calculate maximum profit over all destinations
    double maxProfit = block.getProfit(0);
    for (int d=1; d < prob.getnDestination(); ++d){
      if (block.getProfit(d) > maxProfit)
        maxProfit = block.getProfit(d);
    }
    upitProfit[b] = maxProfit;

    // for each block get predecessors and add arcs
	  for(int p=0; p<block.getNumPred(); ++p){
	    const int pred = block.getPreds()[p];
      upitGraph.addArc(pred,b);
	  }
  }

  //std::cout << "upitGraph has " << upitGraph.getNumNodes() << " nodes\n";
  //std::cout << "upitProfit has " << upitProfit.size() << " elements" << std::endl;
  MaxClosureFactory mcfactory;

  mcfactory.setOption('b');

  upit = mcfactory(upitGraph);

  upit->setProfit(upitProfit);
}

int UpitSolver::solve(){
  return upit->solve();
}


void UpitSolver::getResidualProfit(std::vector<double> &residualProfit){
  upit->getResidualProfit(residualProfit);
}

int UpitSolver::getClosure(std::vector<int> &closure){
  return upit->getClosure(closure);
}

int UpitSolver::getFixed(std::map<int,int> &fixed){
  std::vector<int> closure(nb,0);
  getClosure(closure);

  for (size_t i = 0; i < nb; ++i){
    //std::cout << closure[i] << std::flush;
    if (closure[i] == 0)
      fixed[i] = 0;
  }

  return true;
}

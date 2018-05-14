#include "MaxClosure_PP.h"

MaxClosure_PP::MaxClosure_PP(const Graph &graph, const std::map <int,int> &fixed,
    MaxClosureFactory mcfactory)
                              : MaxClosure_Base(graph) {
  // establish status vector and vector to keep track of in/out degrees
  // status key:  1 - fixed to 1
  //              0 - fixed to 0
  //             -1 - not fixed and not in chain
  //             -2 - not fixed and is in chain
  std::vector<int> status(graph.getNumNodes(), -1);
  std::vector<int> inDegree(status.size(),0);
  //std::cout << "status size: " << status.size() << std::endl;
  std::vector<int> outDegree(status.size(),0);
  _status=0;			// everything OK
  // initialise the degree arrays
  for (size_t i = 0; i < status.size(); i++){
    inDegree[i] = graph.getNode(i)->getInDegree();
    outDegree[i] = graph.getNode(i)->getOutDegree();
  }

  // iterate over fixed map to find nodes to remove
  for (std::map <int,int>::const_iterator it=fixed.begin(); it!=fixed.end(); ++it){

    //std::cout << "key: " << it->first << ", value: " << it->second << std::endl;
    // infeasible solution caused by fixing node
    if (status[it->first] > -1 && status[it->first] != it->second){
	_status = -1;
	return;
    }

    // call recursive function to find nodes affected
    // by fixing (and update in/out degrees)
    // it->second > 0 means go backward and fix all predecessors else forward
    visitDFS_PP(it->first, status, inDegree, outDegree, it->second, it->second);
    if(_status != 0) return;
  }

  // find chains
  std::vector<int> chain; // don't want to re-allocate memory each time
  chain.reserve(10);
  for (int i =0; i < inDegree.size(); i++){
    if (inDegree[i] == 1 && outDegree[i] == 1 && status[i] < 0){
      const Node *curr = graph.getNode(i);
      chain.clear();
      // search backwards in chain
      findChain(curr->getID(), inDegree, outDegree, status, chain, true);
      // remove last value to avoid doubling up
      chain.pop_back();
      // search forwards in chain
      findChain(curr->getID(), inDegree, outDegree, status, chain, false);

      // add current chain to vector of chains
      chains.push_back(chain);
    }
  }

  // build reduced graph
  // Graph reducedGraph = Graph();

  int nodeCount = 0;

  // reserve memory for pp_nodes vector
  pp_nodes.reserve(status.size());

  // if node wasnt in area affected by fixing, create new
  // node and map to problem ID
  for (int i = 0; i < status.size(); i++){
    PP_Node newPP(reducedGraph);
    newPP.status = status[i];
    // if not fixed or in chain, add to reduced graph
    if (status[i] == -1){
      newPP.nodeID = nodeCount++;
    } else {
	    newPP.nodeID = -1; // this gets overriden in setProfit() for chains
    }
    // add pre-processing node to vector
    pp_nodes.push_back(newPP);
  }
  reducedGraph.resizeNodes(nodeCount); // create all nodes in one go

  // update the class member that stores number of nodes in reduced graph
  reducedNumNodes = nodeCount;

  // loop through all nodes and copy appropriate arcs to reduce graph
  for (int i = 0; i < status.size(); i++){
    if (status[i] == -1){
      const Node *curr = graph.getNode(i);
      for (int j = 0; j < curr->getOutDegree(); j++){
        int tgt = curr->getOutArc(j)->getTgtID();
        if (status[tgt] == -1){
          reducedGraph.addArc(pp_nodes[i].nodeID, pp_nodes[tgt].nodeID);
        }
      }
    }
  }

  // update chains in pp_nodes
  for (int i = 0; i < chains.size(); i++){
    // add arc for the chain
      reducedGraph.addArc(pp_nodes[chains[i][0]].nodeID,pp_nodes[chains[i].back()].nodeID);
      for (int j = 1; j < chains[i].size()-1; j++){
	  pp_nodes[chains[i][j]].chain = &chains[i];
      }
  }

  // build reduced problem with reduced graph
  reducedProblem = mcfactory(reducedGraph);
}

MaxClosure_PP::~MaxClosure_PP(){
	if (reducedProblem != NULL) delete reducedProblem;
	reducedProblem = NULL;
}

// run solver for reduced problem
int MaxClosure_PP::solve(){ return reducedProblem->solve(); }

// recursive function to find chains in graph
void MaxClosure_PP::findChain(NodeID curr, std::vector<int> &inDegree,
                              std::vector<int> &outDegree, std::vector<int> &status,
                              std::vector<int> &chain, bool backwards){

  // update status and in/out degrees
  status[curr] = -2;
  inDegree[curr] = 0;
  outDegree[curr] = 0;

  if (backwards){
    // add current ID to beginning of chain vector
/*********************************
    for (int i = 0; i < curr->getInDegree(); ++i){
      int srcID = curr->getInArc(i)->getSrcID();
      // check node "exists" (i.e., not fixed or in a chain)
      if ((outDegree[srcID] > 0 || inDegree[srcID] > 0) && status[srcID] < 0){
        const Node *src = curr->getInArc(i)->getSrc();
        // if previous node is part of chain, recurse
        if (inDegree[src->getID()] == 1 && outDegree[src->getID()] == 1){
          findChain(src, inDegree, outDegree, status, chain, backwards);
        }
        // otherwise add to beginning of chain vector
        else {
		      chain.push_back(src->getID());
        }
        break;
      }
    }
    chain.push_back(curr->getID());
  }
  else{
    // add current ID to end of chain vector
    chain.push_back(curr->getID());
    for (int i = 0; i < curr->getOutDegree(); ++i){
      int tgtID = curr->getOutArc(i)->getTgtID();
      // check node "exists" (i.e., not fixed or in a chain)
      if ((outDegree[tgtID] > 0 || inDegree[tgtID] > 0) && status[tgtID] < 0){
        const Node *tgt = curr->getOutArc(i)->getTgt();
        // if previous node is part of chain, recurse
        if (inDegree[tgt->getID()] == 1 && outDegree[tgt->getID()] == 1){
          findChain(tgt, inDegree, outDegree, status, chain, backwards);
        }
        // otherwise add to end of chain vector
        else {
		      chain.push_back(tgt->getID());
        }
        break;
      }
**********************************/
    //const Node *src = curr->getInArc(0)->getSrc();
      NodeID srcID=_graph->getNode(curr)->getInArc(0)->getSrcID();
    // if previous node is part of chain, recurse - otherwise add to beginning of chain vector
    if (inDegree[srcID] == 1 && outDegree[srcID] == 1){
      findChain(srcID, inDegree, outDegree, status, chain, backwards);
    } else {
		chain.push_back(srcID);
    }
		chain.push_back(curr);
  }else{
    // add current ID to end of chain vector
    chain.push_back(curr);
    NodeID tgtID = _graph->getNode(curr)->getOutArc(0)->getTgtID();
    // if next node is part of chain, recurse - otherwise add to end of chain vector
    if (inDegree[tgtID] == 1 && outDegree[tgtID] == 1){
      findChain(tgtID, inDegree, outDegree, status, chain, backwards);
    }else{
      chain.push_back(tgtID);
    }

//    const Node *tgt = curr->getOutArc(0)->getTgt();
//    // if next node is part of chain, recurse - otherwise add to end of chain vector
//    if (inDegree[tgt->getID()] == 1 && outDegree[tgt->getID()] == 1){
//      findChain(tgt, inDegree, outDegree, status, chain, backwards);
//    }
//    else{
//      //if (status[tgt->getID()] < 0)
//        chain.push_back(tgt->getID());
//    }
  }
}

// function to convert full profit vector to reduced profit
void MaxClosure_PP::setProfit(const std::vector<double> &profit){

  std::vector<double> reducedProfit(reducedNumNodes,0);
  //std::vector<int> status(profit.size(),0);

  // if node exists in the reduced graph, update its profit
  for (int i = 0; i < profit.size(); i++){
    if (pp_nodes[i].status == -1)
      reducedProfit[pp_nodes[i].nodeID] = profit[i];
  }

  // find splitting point in chain and update profits accordingly
  for (int i = 0; i < chains.size(); i++){
    double sum = 0;
    double maxVal = 0,remainVal=0;
    int splitIdx = 0;
    for (int j=1; j < chains[i].size()-1; j++){
      // set node as visited to avoid double handling chains
      //status[chains[i][j]] = 1; // not required
      sum += profit[chains[i][j]];
      if (sum > maxVal){
        maxVal = sum;
        splitIdx = j;
        remainVal=0;
      }else
        remainVal += profit[chains[i][j]];
    }

    // set sum of profits to either side of chain
    const int startNd=pp_nodes[chains[i][0]].nodeID,
		endNd=pp_nodes[chains[i].back()].nodeID;
    reducedProfit[startNd] += maxVal;
    reducedProfit[endNd] += remainVal; //sum-maxVal;

    // assign chain nodes to either side of split
    for(int j = 1; j <= splitIdx; ++j)
	pp_nodes[chains[i][j]].nodeID = startNd;
    for (int j = splitIdx+1; j < chains[i].size()-1; j++)
	pp_nodes[chains[i][j]].nodeID = endNd;
  }

  // call setProfit on reduced problem
  reducedProblem->setProfit(reducedProfit);
}

// function to convert a reduced solution to a full solution
int MaxClosure_PP::getClosure(std::vector<int> &sol){
  // clear input solution
  for(size_t i=0;i<sol.size();++i) sol[i] = 0;

  std::vector<int> reducedSol(reducedNumNodes,0);

  // find reduced closure
  reducedProblem->getClosure(reducedSol);

  // copy values from reduced solution to full solution
  for (size_t i=0; i<sol.size(); i++){
      sol[i] = (pp_nodes[i].status <= -1) ? // read solution from reduced graph
	       reducedSol[pp_nodes[i].nodeID]
	       : pp_nodes[i].status; // in case it's 0/1 for fixed nodes
  }
  return true;
}

void MaxClosure_PP::getResidualProfit(std::vector<double> &residualProfit){
  // clear input solution
  for(size_t i=0;i<residualProfit.size();++i) residualProfit[i] = 0;

  std::vector<double> reducedResidualProfit(reducedNumNodes,0.0);

  // find reduced closure
  reducedProblem->getResidualProfit(reducedResidualProfit);

  // copy values from reduced solution to full solution
  for (size_t i=0; i<residualProfit.size(); i++){
      residualProfit[i] = (pp_nodes[i].status <= -1) ? // read solution from reduced graph
	       reducedResidualProfit[pp_nodes[i].nodeID]
	       : pp_nodes[i].status; // in case it's 0/1 for fixed nodes
  }
}

// Recursive function to perform depth first search
void MaxClosure_PP::visitDFS_PP(int curr, std::vector<int> &status,
                                std::vector<int> &inDegree,
                                std::vector<int> &outDegree,
                                int value, bool backwards){
  // Mark the current node status
  if(status[curr] != -1 && status[curr] != value){
      _status = -1;
      return;
  }
  status[curr] = value;

  // set degree of current node to 0
  inDegree[curr] = 0;
  outDegree[curr] = 0;

  if (backwards){
    for (int i = 0; i < _graph->getNode(curr)->getOutDegree(); i++){
      int next = _graph->getNode(curr)->getOutArc(i)->getTgtID();
      if (inDegree[next] > 0)
        inDegree[next] -= 1; // subtract one from inDegree of next node
    }
    for (int i = 0; i < _graph->getNode(curr)->getInDegree(); i++){
      int prev = _graph->getNode(curr)->getInArc(i)->getSrcID();
      if (status[prev] < 0)
        visitDFS_PP(prev, status, inDegree, outDegree, value, backwards);
    }
  } else {
    for (int i = 0; i < _graph->getNode(curr)->getInDegree(); i++){
      int prev = _graph->getNode(curr)->getInArc(i)->getSrcID();
      if (outDegree[prev] > 0){
        outDegree[prev] -= 1;
      }
    }
    for (int i = 0; i < _graph->getNode(curr)->getOutDegree(); i++){
      int next = _graph->getNode(curr)->getOutArc(i)->getTgtID();
      if (status[next] < 0){
        visitDFS_PP(next, status, inDegree, outDegree, value, backwards);
      }
    }
  }
}

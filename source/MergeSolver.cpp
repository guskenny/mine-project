#include "MergeSolver.h"

MergeSolver::MergeSolver(const SettingsHandler sh, SinglePModel *base_model, BranchNode_info &merged, std::vector<std::vector<int> > &groups, std::vector<int> &group_map) : sh(sh), merged(merged), groups(groups), group_map(group_map) {

  nB = base_model->getNBlock();
  nG = groups.size();
  nR = base_model->getnResources();
  nT = base_model->getNPeriod();
  rate = base_model->getDiscountRate();

  group_profits = std::vector<double>(nG, 0.0);
  group_resources = std::vector<std::vector<double> > (nG, std::vector<double> (nR, 0.0));

  initMergeModel(base_model);

  std::cout << "\nMerge Solver initialised" << std::endl;
}

void MergeSolver::initMergeModel(SinglePModel *base_model){

  // resize the merge graph with number of groups
  merge_graph.resizeNodes(nG);

  // set up matrix to keep track of which groups are connected to what
  std::vector<std::vector<bool> > adj_matrix(nG,std::vector<bool>(nG,false));


  for (int b = 0; b < nB; ++b){
    const Block & block = base_model->getBlock(b);

    // build graph
    const Node* curr = base_model->graph.getNode(b);
    std::vector<Node*> connected = std::vector<Node*>();
    const int num_connected = curr->getConnected(BACKWARD, connected);

    for (int t = 0; t < nT; ++t){
      // iterate over all connected nodes to current node and add group arcs
      for (size_t node = 0; node < num_connected; ++node){
        int node_idx = connected[node]->getID();
        if (group_map[nB*t + b] != group_map[nB*t + node_idx]){
          if (!adj_matrix[group_map[nB*t + node_idx]][group_map[nB*t + b]]){
            merge_graph.addArc(group_map[nB*t + node_idx], group_map[nB*t + b]);
            adj_matrix[group_map[nB*t + node_idx]][group_map[nB*t + b]] = true;
          }
        }
      }

      if (t < nT-1){
        group_profits[group_map[(t*nB)+b]] += (block.getProfit(0)/pow(1+rate,t) - block.getProfit(0)/pow(1+rate,t+1));

        // add arcs to next period groups (arrows point backward in time - if block is 1 at time t, then block must be 1 at time t+1)
        if (group_map[nB*t + b] != group_map[nB*(t+1) + b]){
          if (!adj_matrix[group_map[nB*(t+1) + b]][group_map[nB*t + b]]){
            merge_graph.addArc(group_map[nB*(t+1) + b], group_map[nB*t + b]);
            adj_matrix[group_map[nB*(t+1) + b]][group_map[nB*t + b]] = true;
          }
        }
      }
      else{
        group_profits[group_map[(t*nB)+b]] += (block.getProfit(0)/pow(1+rate,t));
      }
    }
  }
}

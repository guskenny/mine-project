/*************************************************************************** 
                          Implementation of boostMaxFlow class 
                         ------------------------------------------- 
    last modified   : 21/6/2016 
    copyright       : (C) 2016 by Angus Kenny
    libraries		    : . 
    description		  : methods for boostMaxFlow class 
                      (THIS WAS FOR TESTING PURPOSES ONLY!!)
 ***************************************************************************/ 
 
/*************************************************************************** 
 *                                                                         * 
 *   This program is free software; you can redistribute it and/or modify  * 
 *   it under the terms of the GNU General Public License as published by  * 
 *   the Free Software Foundation; either version 2 of the License, or     * 
 *   (at your option) any later version.                                   * 
 *                                                                         * 
 ***************************************************************************/ 

#include "../include/boostMaxFlow.h" 

//using namespace boost;

// constructor for BoostMaxFlow object
BoostMaxFlow::BoostMaxFlow(){
   
  capacity = boost::get(boost::edge_capacity, g);
  rev = boost::get(boost::edge_reverse, g);
  residual_capacity = boost::get(boost::edge_residual_capacity, g);
}

// build boost network from common network class
void BoostMaxFlow::buildNetwork(Network *network){
  long n_nodes = network->getNumNodes();
  long n_arcs = network->getNumArcs();

  // add nodes to graph
  for (long i=0; i < n_nodes; i++)
    verts.push_back(add_vertex(g)); 

  // define source and target nodes
  s = verts[0];
  t = verts[n_nodes-1];

  // add arcs to graph
  for (long i=0; i < n_arcs; i++){
    NetArc *currArc = network->getArc(i);
    int tail = currArc->getSrcID();
    int head = currArc->getTgtID();
    edge_descriptor e1, e2;
    bool in1, in2;
    tie(e1, in1) = add_edge(verts[tail], verts[head], g);
    tie(e2, in2) = add_edge(verts[head], verts[tail], g);

    // add capacities
    capacity[e1] = currArc->getWeight();
    capacity[e2] = 0;
    rev[e1] = e2;
    rev[e2] = e1;
  }
}

// build boost network from common network class
void BoostMaxFlow::buildNetworkNoST(Network *network, std::vector<double> *profits){
  long n_nodes = network->getNumNodes();
  long n_arcs = network->getNumArcs();

  S = n_nodes;
  T = n_nodes+1;

  width = network->getWidth();
  depth = network->getDepth();

  // add nodes to graph
  for (long i=0; i < n_nodes+2; i++)
    verts.push_back(add_vertex(g)); 

  // define source and target nodes
  s = verts[S];
  t = verts[T];

  // add arcs to graph
  for (long i=0; i < n_arcs; i++){
    NetArc *currArc = network->getArc(i);
    int tail = currArc->getSrcID();
    int head = currArc->getTgtID();
    edge_descriptor e1, e2;
    bool in1, in2;
    tie(e1, in1) = add_edge(verts[tail], verts[head], g);
    tie(e2, in2) = add_edge(verts[head], verts[tail], g);

    // add capacities
    capacity[e1] = 1000;
    capacity[e2] = 0;
    rev[e1] = e2;
    rev[e2] = e1;
  }
  
  // add source and sink arcs
  for (int i = 0; i < profits->size(); i++){
    edge_descriptor e1, e2;
    bool in1, in2;
    int tail, head;
    
    if ((*profits)[i] > 0){
      tail = S;
      head = i;
    } 
    else{
      tail = i;
      head = T;
    }

    tie(e1, in1) = add_edge(verts[tail], verts[head], g);
    tie(e2, in2) = add_edge(verts[head], verts[tail], g);

    // add capacities
    capacity[e1] = abs((*profits)[i]);
    capacity[e2] = 0;
    rev[e1] = e2;
    rev[e2] = e1;
  }
}

void BoostMaxFlow::visitDFS(int curr, std::vector<int> &closure)
{
  // Mark the current node as visited and print it
  closure[curr] = 1;
  
  // Recur for all the vertices adjacent to this vertex
  boost::graph_traits<Graph>::out_edge_iterator ei, e_end;
  
  for (boost::tie(ei, e_end) = out_edges(verts[curr], g); ei != e_end; ++ei)
    if (residual_capacity[*ei] > 0){
      int next = target(*ei, g);
      if (!closure[next])
        visitDFS(next,closure); 
    }
}

void BoostMaxFlow::printFlowOutput(long flow){
  
  std::cout << "c  The total flow:" << std::endl;
  std::cout << "s " << flow << std::endl << std::endl;

  std::cout << "c flow values:" << std::endl;
 
  boost::graph_traits<Graph>::vertex_iterator u_iter, u_end;
  boost::graph_traits<Graph>::out_edge_iterator ei, e_end;
  
  for (boost::tie(u_iter, u_end) = vertices(this->g); u_iter != u_end; ++u_iter)
    for (boost::tie(ei, e_end) = out_edges(*u_iter, this->g); ei != e_end; ++ei)
      if (this->capacity[*ei] > 0 ){ 
        std::cout << "f " << *u_iter << " " << target(*ei, this->g) << " " 
                  << this->residual_capacity[*ei] << std::endl;
     }
}

void BoostMaxFlow::printClosure(std::vector<int> &closure){

  std::cout <<"\nCLOSURE VECTOR:\n";

  for (int i=0; i < closure.size(); i++){
    std::cout << closure[i] << " ";
  }

  std::cout << "\n\nMAX CLOSURE NODES:\n";
  // print out max closure nodes
  for (int i=0; i < closure.size(); i++){
    if (closure[i] == 1)
      std::cout << i << " ";
  }
  std::cout << std::endl;
}

void BoostMaxFlow::printClosureGraph(std::vector<int> &closure, int depth, int width){
  char green[] = { 0x1b, '[', '1', ';', '3', '2', 'm', 0 };
  char red[] = { 0x1b, '[', '1', ';', '3', '1', 'm', 0 };
  char blue[] = { 0x1b, '[', '1', ';', '3', '4', 'm', 0 };
  char normal[] = { 0x1b, '[', '0', ';', '3', '9', 'm', 0 };

  std::cout << "\nCLOSURE GRAPH\n\n";

  int i = 0;
  for (int j=0; j < depth; j++){
    for (int k=0; k < width; k++){
      if (closure[i] == 1)
        std::cout << blue << closure[i]+1 << " ";
      else if (closure[i] == 0)
        std::cout << red << closure[i]+1 << " ";
      else if (closure[i] == -1)
        std::cout << green << closure[i]+1 << " ";
      i++;
    }
    std::cout << std::endl;
  }
  std::cout << normal << "\n*** All values have had 1 added to them for display ***" << std::endl;
}
bool BoostMaxFlow::getClosure(std::vector<int> &sol){ 
	visitDFS(S,sol); 
	return true; 
}
void BoostMaxFlow::pushRelabel(){ 
	boost::push_relabel_max_flow(g, s, t); 
}

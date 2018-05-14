/*************************************************************************** 
                            A Class for Boost Push Relabel 
                         ---------------------------------- 
    last modified   : 21/6/2016 
    copyright       : (C) 2016 by Angus Kenny 
    libraries		    : . 
    description		  : contains data structure for boost algorithms
                      (USED FOR TESTING PURPOSES ONLY!!!)
 ***************************************************************************/ 
 
/*************************************************************************** 
 *                                                                         * 
 *   This program is free software; you can redistribute it and/or modify  * 
 *   it under the terms of the GNU General Public License as published by  * 
 *   the Free Software Foundation; either version 2 of the License, or     * 
 *   (at your option) any later version.                                   * 
 *                                                                         * 
 ***************************************************************************/ 

#include <boost/config.hpp>
#include <iostream>
#include <string>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/read_dimacs.hpp>
#include <boost/graph/graph_utility.hpp>

#include "network.h"

class BoostMaxFlow{
  private:
    // define all the boost types  
    typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> Traits;
    typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, 
      boost::property<boost::vertex_name_t, std::string>,
      boost::property<boost::edge_capacity_t, long,
        boost::property<boost::edge_residual_capacity_t, long,
          boost::property<boost::edge_reverse_t, Traits::edge_descriptor> > >
    > Graph;
  
    typedef typename boost::graph_traits<Graph>::vertices_size_type vertices_size_type;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
  
    Graph g;

    boost::property_map<Graph, boost::edge_capacity_t>::type capacity;
    boost::property_map<Graph, boost::edge_reverse_t>::type rev;
    boost::property_map<Graph, boost::edge_residual_capacity_t>::type residual_capacity;
  
    Traits::vertex_descriptor s, t;
  
    std::vector<vertex_descriptor> verts;

    int width;
    int depth;

    int S;
    int T;

  public:
    BoostMaxFlow();
    ~BoostMaxFlow(){};
    void buildGraph(Network *network, std::vector<double> *profits);
    void buildNetwork(Network *network);
    void buildNetworkNoST(Network *network, std::vector<double> *profits);
    void pushRelabel();
    bool getClosure(std::vector<int> &sol); 
    void visitDFS(int curr, std::vector<int> &closure);
    void printFlowOutput(long flow);
    void printClosure(std::vector<int> &closure);
    static void printClosureGraph(std::vector<int> &closure, int depth, int width);
};

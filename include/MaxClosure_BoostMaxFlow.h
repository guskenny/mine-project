/*************************************************************************** 
                     A template class for boost max flow 
                    ------------------------------------- 
    last modified   : 21/6/2016 
    copyright       : (C) 2016 by Angus Kenny 
    libraries		    : . 
    description		  : contains data structure for boost algorithms, inherits
                      from MaxClosure_Base class. Class is a template class
                      which allows different boost algorithms to be used.
                      
                      Use the classes:
                      
                      MaxClosure_BoostMaxFlow_PL(graph) for push-relabel
                      MaxClosure_BoostMaxFlow_BK(graph) for Boykov-Kolmogorov
                      MaxClosure_BoostMaxFlow_EK(graph) for Edmonds-Karp

                      Function definitions have not been separated from the
                      header file because it is a template.

 ***************************************************************************/ 
 
/*************************************************************************** 
 *                                                                         * 
 *   This program is free software; you can redistribute it and/or modify  * 
 *   it under the terms of the GNU General Public License as published by  * 
 *   the Free Software Foundation; either version 2 of the License, or     * 
 *   (at your option) any later version.                                   * 
 *                                                                         * 
 ***************************************************************************/ 

#ifndef MaxClosure_BoostMaxFlow_H
#define MaxClosure_BoostMaxFlow_H

#include "boostGraphTypes.h"
#include "MaxClosure_Base.h"
   
template <class GraphType>
class MaxClosure_BoostMaxFlow : public MaxClosure_Base {
  public: 
    // Boost graph type, defined by specific derived class
    typedef GraphType BoostGraph;
  
  protected:
    // Boost type definitions
    typedef typename boost::graph_traits<BoostGraph>::vertices_size_type 
      vertices_size_type;
    typedef typename boost::graph_traits<BoostGraph>::vertex_descriptor 
      vertex_descriptor;
    typedef typename boost::graph_traits<BoostGraph>::edge_descriptor 
      edge_descriptor;

    typename boost::property_map < BoostGraph, boost::edge_capacity_t >::type capacity;
    typename boost::property_map < BoostGraph, boost::edge_residual_capacity_t >::type residual_capacity;
    typename boost::property_map < BoostGraph, boost::edge_reverse_t >::type rev;
    
    // constant "dummy infinity" value for uncapacitated arcs
    static constexpr double DUMMY_INF = 1.0e18;
    BoostGraph g;
    Traits::vertex_descriptor s, t;
    std::vector<vertex_descriptor> verts;
    long n_nodes; // number of nodes in problem instance (excluding S,T)
    long n_arcs; // number of arcs in problem instance (excluding S,T arcs)
    long S;
    long T;

  public:
    // Constructor
    MaxClosure_BoostMaxFlow(const Graph &graph) : MaxClosure_Base(graph) {
       
      capacity = get(boost::edge_capacity, g);
      rev = get(boost::edge_reverse, g);
      residual_capacity = get(boost::edge_residual_capacity, g);
    
      n_nodes = graph.getNumNodes();
      n_arcs = graph.getNumArcs();
    
      S = n_nodes;
      T = n_nodes+1;
    
      // add nodes to graph
      for (long i=0; i < n_nodes+2; i++)
        verts.push_back(boost::add_vertex(g)); 
    
      // define source and target nodes
      s = verts[S];
      t = verts[T];
    
      // add arcs to graph
      for (long i=0; i < n_arcs; i++){
        const Arc *currArc = graph.getArc(i);
        int tail = currArc->getSrcID();
        int head = currArc->getTgtID();
        edge_descriptor e1, e2;
        bool in1, in2;
        tie(e1, in1) = add_edge(verts[head], verts[tail], g);
        tie(e2, in2) = add_edge(verts[tail], verts[head], g);
    
        // add capacities
        capacity[e1] = DUMMY_INF;
        capacity[e2] = 0;
        rev[e1] = e2;
        rev[e2] = e1;
      }
	};
	~MaxClosure_BoostMaxFlow(){};
    
    // Function to set arc capacities as per provided profit vector
    void setProfit(const std::vector<double> &profit){
      
      boost::clear_vertex(verts[S],g);
      boost::clear_vertex(verts[T],g);

      // add source and sink arcs
      for (int i = 0; i < profit.size(); i++){
        edge_descriptor e1, e2;
        bool in1, in2;
        int tail, head;
        
        if (profit[i] > 0){
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
        capacity[e1] = abs(profit[i]);
        capacity[e2] = 0;
        rev[e1] = e2;
        rev[e2] = e1;
      }
	};

    // function to compute residual profits after max closure has been solved
    void getResidualProfit(std::vector<double> &residualProfit){
      // make sure residual profits are all zero
      for (size_t i=0;i<residualProfit.size();++i) 
        residualProfit[i] = 0;

      typename boost::graph_traits<BoostGraph>::out_edge_iterator ei, e_end;

      // iterate over all edges coming from the source and store any positive
      // values as blocks with residual profit
      for (boost::tie(ei, e_end) = out_edges(verts[S], g); ei != e_end; ++ei)
        if (residual_capacity[*ei] > 0){
          const int tgt = target(*ei,g);
          residualProfit[tgt] = residual_capacity[*ei];
        }
    };

    // Virtual function, to be implemented by specific class
    virtual int solve()=0;

    // Function to use DFS to get closure from a given solution
    int getClosure(std::vector<int> &sol){
      // add dummy nodes to sol vector
      for(size_t i=0;i<sol.size();++i) sol[i] = 0; // clear solution
      sol.push_back(0);
      sol.push_back(0);
    
      visitDFS(S,sol); 
    
      // remove dummy nodes from sol vector
      sol.pop_back();
      sol.pop_back();
    
      return true;
	};

    // Recursive function to perform depth first search
    void visitDFS(int curr, std::vector<int> &closure){
      // Mark the current node as visited and print it
      std::vector<int> stack;
      stack.reserve(closure.size());
      stack.push_back(curr);
      // Recurse for all the vertices adjacent to this vertex
      typename boost::graph_traits<BoostGraph>::out_edge_iterator ei, e_end;
      closure[curr] = 1;
      while( ! stack.empty() ){
		  curr =  stack.back(); stack.pop_back();
		  for (boost::tie(ei, e_end) = out_edges(verts[curr], g); ei != e_end; ++ei)
			  if (residual_capacity[*ei] > 0){
				  const int next = target(*ei, g);
				  if (!closure[next]){
					  closure[next] = 1;
					  stack.push_back(next);
				  }
			  }
	  }
   }; 


    // get dual solution information:
    //bool getInFlowSolution(std::vector< std::vector<double> > &inFlow){};  
    //inFlow[i][j] is flow from node pred[i][j] to i
    // outFlow is really the same information as inFlow
    // virtual bool getOutFlowSolution(std::vector< std::vector<double> > &outFlow)=0; 
    //outFlow[i][j] is flow from node i to succ[i][j] to i
    //bool getRedCost(std::vector<double> &sol){}; 
    // reduced cost (dual) associated with 0<= x_i <= 1 constraints
   
    
};
#endif


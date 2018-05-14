

/*************************************************************************** 
                          Implementation of LR_Graph class 
                         ------------------------------- 
    last modified   : 21/6/2016 
    copyright       : (C) 2016 by Angus Kenny
                      (based on network.cpp by Dhananjay Thiruvady)
    libraries    : . 
    description  : the graph data structure
***************************************************************************/ 
 
/*************************************************************************** 
 *                                                                         * 
 *   This program is free software; you can redistribute it and/or modify  * 
 *   it under the terms of the GNU General Public License as published by  * 
 *   the Free Software Foundation; either version 2 of the License, or     * 
 *   (at your option) any later version.                                   * 
 *                                                                         * 
 ***************************************************************************/ 

using namespace std; 

#include "lr_graph.h"

/*
  LR_Graph related function implementations
*/


/*LR_Graph::LR_Graph(int nVertices,const vector<vector<int> > &_pred){
  this->pred = _pred;
  }*/


// Method to add node to graph
/*void LR_Graph::addExtNode(int nodeID){
  ExtNode *newExtNode = new ExtNode(nodeID);
  nodes.push_back(newExtNode);
}

// Method to add arc to graph
void LR_Graph::addExtArc(int arcID, int src, int tgt) {
  ExtArc *newExtArc = new ExtArc(arcID, nodes[src], nodes[tgt]);
  arcs.push_back(newExtArc);
  nodes[src]->insertOutExtArc(newExtArc);
  nodes[tgt]->insertInExtArc(newExtArc);
}

// Method to get specific arc in graph (returns NULL pointer if no arc exists)
ExtArc *LR_Graph::getExtArcPair(int srcID, int tgtID) {
  for (int i=0; i < nodes[srcID]->getOutDegree(); i++) 
    if (nodes[srcID]->getOutExtArc(i)->getTgtID() == tgtID)
      return nodes[srcID]->getOutExtArc(i);
  
  return NULL;
}

const vector<int> & LR_Graph::getPreds(int vertex){
  return pred[vertex];
}

// The graph destructor
LR_Graph::~LR_Graph(){
  for(int i=0;i<nodes.size();i++){
    delete nodes[i];
  }
  for(int i=0;i<arcs.size();i++){
    delete arcs[i];
  }
  }*/

/*
  ExtNode related function implementations
*/

// ExtNode class constructor
ExtNode::ExtNode(int id){
  this->id = id;
}

/*
  ExtArc related function implementations
*/

// ExtArc class constructor
ExtArc::ExtArc(int id, ExtNode *src, ExtNode *tgt){
  this->id = id;
  this->src = src;
  this->tgt = tgt;
}

/*
  Net ExtNode related function implementations
*/

// ExtNode class constructor
NetExtNode::NetExtNode(int id, double supply, double value) : ExtNode(id){
  this->id = id;
  this->supply = supply;
  this->value = value;
}


/*
  ExtArc related function implementations
*/

// ExtArc class constructor
NetExtArc::NetExtArc(int id, NetExtNode *src, NetExtNode *tgt, double lb, double ub, double cost){
  this->id;
  this->src = src;
  this->tgt = tgt;
  this->lb = lb;
  this->ub = ub;
  this->cost = cost;
}


// The arc destructor
NetExtArc::~NetExtArc(){
  if(src) delete src;
  if(tgt) delete tgt;
}


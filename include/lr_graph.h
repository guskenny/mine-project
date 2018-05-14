

/*************************************************************************** 
                            A Graph Class 
                         ------------------- 
    last modified   : 21/6/2016 
    copyright       : (C) 2016 by Angus Kenny 
                      (based on network.h by Dhananjay Thiruvady)
    libraries  : . 
    description: contains a data structure for a graph
***************************************************************************/ 
 
/*************************************************************************** 
 *                                                                         * 
 *   This program is free software; you can redistribute it and/or modify  * 
 *   it under the terms of the GNU General Public License as published by  * 
 *   the Free Software Foundation; either version 2 of the License, or     * 
 *   (at your option) any later version.                                   * 
 *                                                                         * 
 ***************************************************************************/ 

#ifndef LR_Graph_H
#define LR_Graph_H

#include <iostream>
#include <new>
#include <vector>
#include <string>
#include <set>
#include <cstdio>
#include <cstring>

using namespace std;

class ExtArc;
class ExtNode;

class ExtNode{
 private:
  int id;
  vector<ExtArc*> inExtArcs;
  vector<ExtArc*> outExtArcs;

 public:
  ExtNode(int id);
  void insertInExtArc(ExtArc *in){ inExtArcs.push_back(in);};
  void insertOutExtArc(ExtArc *out){ outExtArcs.push_back(out);};
  int getID() {  return id;};
  int getInDegree() {return inExtArcs.size();};
  int getOutDegree() {return outExtArcs.size();};
  int getDegree() {return this->getInDegree() + this->getOutDegree();};
  ExtArc *getInExtArc(int index){ return inExtArcs[index];};
  ExtArc *getOutExtArc(int index){ return outExtArcs[index];};
  ~ExtNode(){};
};

class ExtArc{
 protected:
  int id;
  ExtNode *src;
  ExtNode *tgt;
 public:
  ExtArc(int id, ExtNode *src, ExtNode *tgt);
  ExtNode *getSrc(){ return src;};
  ExtNode *getTgt(){ return tgt;};
  int getSrcID() const { return src->getID();};
  int getTgtID()const { return tgt->getID();};
  int getID(){ return id;};
  ~ExtArc(){};
};

/*class LR_Graph{
 protected:
  vector<ExtNode*> nodes;
  vector<ExtArc*> arcs;
  vector<vector<int> > pred;

 public:
  LR_Graph(int nVertices,const vector<vector<int> > &_pred);
  ExtNode *getExtNode(int nodeID){ return nodes[nodeID];};
  ExtArc *getExtArc(int arcID) const {return arcs[arcID];};
  ExtArc *getExtArcPair(int srcID, int tgtID);
  long getNumExtNodes() const { return nodes.size();};
  long getNumExtArcs() const { return arcs.size();};
  void addExtNode(int nodeID); 
  void addExtArc(int arcID, int src, int tgt);
  const vector<int> & getPreds(int vertex);
  ~LR_Graph();
  };*/


/*
  Network related data structures
*/

class NetExtArc;
class NetExtNode;

// A simple node class
// Stores values of the output of the MIP in val

class NetExtNode : public ExtNode{
 private:
  int id;
  double supply;
  double value;
 public:
  NetExtNode(int id, double supply = 0.0, double value = 0.0);
  void setSupply(double supply){  this->supply = supply;};
  void setValue(double val){    this->value = val;};
  double getSupply(){   return supply;};
  int getID() {  return id;};
  double getValue(){  return value;};
  ~NetExtNode(){};
};

// A simple ExtArc class
// Stores values of the output of the MIP in val

class NetExtArc{
 private:
  int id;
  NetExtNode *src;
  NetExtNode *tgt;
  double lb;
  double ub;
  double cost;
  double value;
 public:
  NetExtArc(int id, NetExtNode *src, NetExtNode *to, double lb, double ub, double cost);
  NetExtNode *getSrc(){  return src;};
  NetExtNode *getTgt(){  return tgt;};
  void setCost(double cost){  this->cost = cost;};
  double getLb(){ return lb;};
  double getUb(){ return ub;};
  double getCost(){  return cost;};
  int getID(){  return id;};
  ~NetExtArc();
};

#endif


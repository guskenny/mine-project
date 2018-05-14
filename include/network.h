/*************************************************************************** 
                            A Network Class 
                         ------------------- 
    last modified   : 21/6/2016 
    copyright       : (C) 2016 by Angus Kenny 
                      (based on work by Dhananjay Thiruvady)
    libraries		: . 
    description		: contains a data structure for a network
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

#ifndef Network_H
#define Network_H

#include <iostream>
#include <new>
#include <vector>
#include <string>
#include <set>
#include <cstdio>
#include <cstring>

#include "graph.h"

class NetArc;
class NetNode;

class Network{
 public:
  Network();
  ~Network();
  
  NetNode *getNode(int nodeID);
  NetArc *getArc(int arcID);
  long getNumNodes();
  long getNumArcs();
  int getWidth();
  int getDepth();
  void addNode(int nodeID); 
  void addArc(int arcID, int src, int tgt, double weight);
  void DFSUtil(NetNode *curr, bool visited[]);
  void readDimacs();
  std::vector<double> readDimacsNoST();
  void printNet();
  std::vector<NetNode*> getNodes();
  std::vector<NetArc*> getArcs();
  Graph makeGraph();

private:
  std::vector<NetNode*> nodes;
  std::vector<NetArc*> arcs;
  int width;
  int depth;
 
};

class NetNode{
 public:
  NetNode(int id);
  ~NetNode();
  
  void insertInArc(NetArc *in);
  void insertOutArc(NetArc *out);
  void setClosure(bool closureStatus);
  bool inClosure();
  int getID();
  int getInDegree();
  int getOutDegree();
  int getDegree();
  NetArc *getInArc(int index);
  NetArc *getOutArc(int index);


 private:
  int id;
  bool isClosure;
  std::vector<NetArc*> inArc;
  std::vector<NetArc*> outArc;

};

class NetArc{
 public:
  NetArc(int id, NetNode *src, NetNode *tgt, double weight);
  ~NetArc();
  NetNode *getSrc();
  NetNode *getTgt();
  int getSrcID();
  int getTgtID();
  void setWeight(double weight);
  double getWeight();
  int getID();

private:
  int id;
  NetNode *src;
  NetNode *tgt;
  double weight;
};

#endif

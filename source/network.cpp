/***************************************************************************
                          Implementation of Network class
                         -------------------------------------------
    last modified   : 21/6/2016
    copyright       : (C) 2016 by Angus Kenny
                      (based on work by Dhananjay Thiruvady)
    libraries		    : .
    description		  : the network data structure
                      (THIS WAS FOR TESTING PURPOSES ONLY!!!)
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

/*
  Network related function implementations
*/
Network::Network(){}
// The network destructor
Network::~Network(){
  for(int i=0;i<nodes.size();i++){
    delete nodes[i];
  }
  for(int i=0;i<arcs.size();i++){
    delete arcs[i];
  }
}

NetNode *Network::getNode(int nodeID){ return nodes[nodeID];}
NetArc *Network::getArc(int arcID) {return arcs[arcID];}
long Network::getNumNodes(){ return nodes.size();}
long Network::getNumArcs() { return arcs.size();}
int Network::getWidth() { return width;}
int Network::getDepth() { return depth;}

void Network::addNode(int nodeID){
  NetNode *newNode = new NetNode(nodeID);
  nodes.push_back(newNode);
}

void Network::addArc(int arcID, int src, int tgt, double weight){
  NetArc *newArc = new NetArc(arcID, nodes[src], nodes[tgt], weight);
  arcs.push_back(newArc);
  nodes[src]->insertOutArc(newArc);
  nodes[tgt]->insertInArc(newArc);
}

void Network::DFSUtil(NetNode *curr, bool visited[])
{
  // Mark the current node as visited and print it
  visited[curr->getID()] = true;
//  cout << curr->getID() << " ";
  // Recur for all the vertices adjacent to this vertex
  for (int i = 0; i < curr->getOutDegree(); ++i)
    if (!visited[curr->getOutArc(i)->getTgtID()])
      if (curr->getOutArc(i)->getWeight() > 0)
        DFSUtil(curr->getOutArc(i)->getTgt(), visited);
}

/*
   Method to read a DIMACS file to network object for testing
*/
void Network::readDimacs(){

  char pr_type[3];
  long n_nodes, n_arcs, i, head, tail, cap;
  long arcID = 0;
  std::istream& in = std::cin;
  std::string in_line;

  while (std::getline(in, in_line)) {
    switch (in_line[0]) {
      case 'c':       // skip comment lines
        break;
      case '\n':      // skip empty lines
        break;
      case '\0':      // skip empty lines at end of file
        break;
      case 'n':       // skip node description (info from arcs)
        break;
      case 'p':       // read problem description
        sscanf (in_line.c_str(), "%*c %3s %ld %ld", pr_type,&n_nodes,&n_arcs);
        //std::cout << "pr_type=" << pr_type
        //          << ", n_nodes=" << n_nodes
        //          << ", n_arcs=" << n_arcs
        //          << std::endl;

        for (i=0; i<n_nodes; i++){
          this->addNode(i);
          //std::cout << "node " << i << " added\n";
        }
        break;
      case 'a':       // read arc description
        sscanf ( in_line.c_str(),"%*c %ld %ld %ld",
                               &tail, &head, &cap );
        this->addArc(arcID++, tail-1, head-1, double(cap));
        //std::cout << "arc (" << tail-1 << ", " << head-1 << ") added\n";
        break;
    }
  }
}

std::vector<double> Network::readDimacsNoST(){

  char pr_type[3];
  long width, depth, n_nodes, n_arcs, node, profit, i, head, tail, cap;
  long arcID = 0;
  std::istream& in = std::cin;
  std::string in_line;
  bool dflag = false;

  std::vector<double> profits;

  while (std::getline(in, in_line)) {
    switch (in_line[0]) {
      case 'c':       // skip comment lines
        break;
      case '\n':      // skip empty lines
        break;
      case '\0':      // skip empty lines at end of file
        break;
      case 'n':       // skip node description (info from arcs)
        break;
      case 'p':       // read problem description
        sscanf (in_line.c_str(), "%*c %3s %ld %ld", pr_type,&n_nodes,&n_arcs);
//        std::cout << "pr_type=" << pr_type
//                  << ", n_nodes=" << n_nodes
//                  << ", n_arcs=" << n_arcs
//                  << std::endl;

        for (i=0; i<n_nodes; i++){
//          std::cout << "adding node " << i << "... ";
          this->addNode(i);
          profits.push_back(0);
//          std::cout << "done!\n";
        }
        break;
      case 'd':
        dflag = true;
        sscanf ( in_line.c_str(),"%*c %ld %ld",
            &width, &depth );
//        std::cout << "adding width: " << width << " and depth: " << depth << "... ";
        this->width = width;
        this->depth = depth;
//        std::cout<< "done!\n";
        break;
      case 's':
        sscanf ( in_line.c_str(),"%*c %ld %ld",
            &node, &profit );
//        std::cout << "adding profit " << profit << " to node " << node << "... ";
        profits[node] = profit;
//        std::cout << "done!\n";
        break;
      case 'a':       // read arc description
        sscanf ( in_line.c_str(),"%*c %ld %ld %ld",
            &tail, &head, &cap );
//        std::cout << "adding arc (" << tail << ", " << head << ")... ";
        this->addArc(arcID++, tail, head, double(cap));
//        std::cout << "done!\n";
        break;
    }
  }
  if (!dflag){
    this->width = 0;
    this->depth = 0;
  }

  return profits;
}

void Network::printNet(){
  std::cout << "Printing nodes:\n";
  for(int i=0; i<nodes.size();i++){
    std::cout << nodes[i]->getID() << ": in="
      << nodes[i]->getInDegree() << ", out="
      << nodes[i]->getOutDegree() << ", total="
      << nodes[i]->getDegree() << std::endl;
  }
  std::cout << "\nPrinting arcs:\n";
  for(int i=0; i < arcs.size(); i++){
    std::cout << "(" << arcs[i]->getSrc()->getID()
      << ", " << arcs[i]->getTgt()->getID() << ") - weight="
      << arcs[i]->getWeight() << std::endl;

  }
}
std::vector<NetNode*> Network::getNodes() {return nodes;}
std::vector<NetArc*> Network::getArcs() {return arcs;}

Graph Network::makeGraph(){
  Graph myGraph = Graph();

  for (int i=0; i < nodes.size();i++)
    myGraph.addNode(nodes[i]->getID());

  for (int i=0; i< arcs.size();i++)
    myGraph.addArc(arcs[i]->getID(), arcs[i]->getSrcID(), arcs[i]->getTgtID());

  return myGraph;
}



/*
  Node related function implementations
*/

NetNode::NetNode(int id){
  this->id = id;
  this->isClosure = false;
}
NetNode::~NetNode(){}

void NetNode::insertInArc(NetArc *in){  inArc.push_back(in); }
void NetNode::insertOutArc(NetArc *out){   outArc.push_back(out);}

/*
  Needs to implement a more complex csetting of the cost
  This will come shortly
*/
void NetNode::setClosure(bool closureStatus){
  this->isClosure = closureStatus;
}

bool NetNode::inClosure(){ return isClosure;}
int NetNode::getID() {  return id;}
int NetNode::getInDegree() {return inArc.size();}
int NetNode::getOutDegree() {return outArc.size();}
int NetNode::getDegree() {return this->getInDegree() + this->getOutDegree();}
NetArc *NetNode::getInArc(int index){ return inArc[index];}
NetArc *NetNode::getOutArc(int index){ return outArc[index];}

/*
  Arc related function implementations
*/

NetArc::NetArc(int id, NetNode *src, NetNode *tgt, double weight){
  this->id = id;
  this->src = src;
  this->tgt = tgt;
  this->weight = weight;
}
NetArc::~NetArc(){}
NetNode *NetArc::getSrc(){  return src;}
NetNode *NetArc::getTgt(){  return tgt;}
int NetArc::getSrcID(){  return src->getID();}
int NetArc::getTgtID(){  return tgt->getID();}

/*
  Needs to implement a more complex csetting of the cost
  This will come shortly
*/

void NetArc::setWeight(double weight){
  this->weight = weight;
}
double NetArc::getWeight(){  return weight;}
int NetArc::getID(){  return id;}

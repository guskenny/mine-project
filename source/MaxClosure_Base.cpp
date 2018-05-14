/***************************************************************************
                            MaxClosure_Base.cpp
                         ----------------------------------
    last modified   : 13/5/2008
    copyright       : (C) 2013 by Dhananjay Thiruvady
    libraries	    : .
    description	    : Implementation of functions in MaxClosure_Base
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 
#include <iostream>
#include <fstream>

using namespace std;

#include "../include/MaxClosure_Base.h"
MaxClosure_Base::MaxClosure_Base(const Graph &graph) : _graph(&graph) {}
MaxClosure_Base::~MaxClosure_Base(){}
int MaxClosure_Base::getStatus(){return _status;}

double MaxClosure_Base::calcProfit(const std::vector<double> &profit,const std::vector<int> &soln) const
{ 
	double sum=0; 
	for(size_t i=0;i<soln.size();++i) 
		if(soln[i]) sum += profit[i]; 
	return sum; 
}
bool MaxClosure_Base::isClosure(const std::vector<int> &soln) const
{
  for(int a=0;a<_graph->getNumArcs();++a){
    const Arc *arc=_graph->getArc(a);
    if( soln[arc->getTgtID()] && ! soln[arc->getSrcID()] ){
      std::cerr << "ERROR: invalid solution "
		<< arc->getSrcID() << " -> " << arc->getTgtID()
		<< std::endl;
      return false;
    }
  }
  return true;
}

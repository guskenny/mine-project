/*************************************************************************** 
              A Class for Boost Max Flow (Edmonds-Karp algorithm)
            ------------------------------------------------------- 
    last modified   : 29/6/2016 
    copyright       : (C) 2016 by Angus Kenny 
    libraries		    : . 
    description		  : contains data structure for boost Edmonds-Karp
                      algorithm, inherits from MaxClosure_BoostMaxFlow
                      template class. Class defines BoostGraph type and
                      implements solve()
 ***************************************************************************/ 
 
/*************************************************************************** 
 *                                                                         * 
 *   This program is free software; you can redistribute it and/or modify  * 
 *   it under the terms of the GNU General Public License as published by  * 
 *   the Free Software Foundation; either version 2 of the License, or     * 
 *   (at your option) any later version.                                   * 
 *                                                                         * 
 ***************************************************************************/ 

#ifndef MaxClosure_BoostMaxFlow_EK_H
#define MaxClosure_BoostMaxFlow_EK_H

#include <boost/graph/edmonds_karp_max_flow.hpp>
#include "MaxClosure_BoostMaxFlow.h"
class MaxClosure_BoostMaxFlow_EK : public MaxClosure_BoostMaxFlow <PREKGraph> {
  public:
    MaxClosure_BoostMaxFlow_EK(const Graph &graph);
    ~MaxClosure_BoostMaxFlow_EK();
    int solve();
};   
#endif


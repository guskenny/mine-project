/*************************************************************************** 
           Type definitions for different boost max flow algorithms
          ----------------------------------------------------------
    last modified   : 29/6/2016 
    copyright       : (C) 2016 by Angus Kenny 
    libraries		    : . 
    description		  : contains different graph types for use with the
                      MaxClosure_BoostMaxFlow template
 ***************************************************************************/ 
 
/*************************************************************************** 
 *                                                                         * 
 *   This program is free software; you can redistribute it and/or modify  * 
 *   it under the terms of the GNU General Public License as published by  * 
 *   the Free Software Foundation; either version 2 of the License, or     * 
 *   (at your option) any later version.                                   * 
 *                                                                         * 
 ***************************************************************************/ 

#ifndef BoostGraphTypes_H
#define BoostGraphTypes_H
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/read_dimacs.hpp>
#include <boost/graph/graph_utility.hpp>

// Common "Traits" type
typedef boost::adjacency_list_traits < boost::vecS, 
        boost::vecS, boost::directedS > Traits;

// Graph type for boost Boykov-Kolmogorov max flow algorithm
typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::directedS,
// boost::property < boost::vertex_name_t, std::string, //AE: don't need name
        boost::property < boost::vertex_index_t, int,
        boost::property < boost::vertex_color_t, boost::default_color_type,
        boost::property < boost::vertex_distance_t, long,
        boost::property < boost::vertex_predecessor_t, Traits::edge_descriptor
			  > > > >, // >,
        boost::property < boost::edge_capacity_t, long,
        boost::property < boost::edge_residual_capacity_t, long,
        boost::property < boost::edge_reverse_t, 
        Traits::edge_descriptor > > >
			  ,boost::vecS> BKGraph;

// Graph type for boost push relabel and Edmonds-Karp max flow algorithm
typedef boost::adjacency_list<boost::listS, boost::vecS,        
			      boost::directedS,
			      boost::property<boost::vertex_name_t, std::string>, boost::property<boost::edge_capacity_t, 
        long, boost::property<boost::edge_residual_capacity_t, 
        long, boost::property<boost::edge_reverse_t, 
        Traits::edge_descriptor> > > > PREKGraph;

#endif

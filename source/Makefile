########################################################################### 
#              Boost push relabel algorithm to solve max closure
#              ------------------------------------------------- 
#    last modified   : 21/6/2016 
#    copyright       : (C) 2016 by Angus Kenny 
#    libraries		    : boost
#    description		  : uses boost push relabel algorithm to 
#                       solve max closure
# 
#                                                                          
#   This program is free software; you can redistribute it and/or modify   
#   it under the terms of the GNU General Public License as published by   
#   the Free Software Foundation; either version 2 of the License, or      
#   (at your option) any later version.                                    
#                                                                          
########################################################################### 


CC=g++
CFLAGS=-ansi -c -std=c++11

all: max_closure gen_graph

gen_graph: gen_graph.o
	$(CC) gen_graph.o -o gen_graph

gen_graph.o: gen_graph.cpp
	$(CC) $(CFLAGS) gen_graph.cpp

#max_closure: max_closure_no_st.o boostMaxFlow.o network.o graph.o MaxClosure_BoostMaxFlow.o 
#	$(CC) max_closure_no_st.o boostMaxFlow.o network.o graph.o MaxClosure_BoostMaxFlow.o -o max_closure

max_closure: max_closure_no_st.o boostMaxFlow.o network.o graph.o MaxClosure_Base.o
	$(CC) max_closure_no_st.o boostMaxFlow.o network.o graph.o MaxClosure_Base.o -o max_closure

max_closure_no_st.o: max_closure_no_st.cpp
	$(CC) $(CFLAGS) max_closure_no_st.cpp

boostMaxFlow.o: boostMaxFlow.cpp
	$(CC) $(CFLAGS) boostMaxFlow.cpp

network.o: network.cpp
	$(CC) $(CFLAGS) network.cpp

graph.o: graph.cpp
	$(CC) $(CFLAGS) graph.cpp

MaxClosure_Base.o: MaxClosure_Base.cpp
	$(CC) $(CFLAGS) MaxClosure_Base.cpp

#MaxClosure_BoostMaxFlow.o: MaxClosure_BoostMaxFlow.cpp
#	$(CC) $(CFLAGS) MaxClosure_BoostMaxFlow.cpp

clean:
	rm -f *.o max_closure gen_graph

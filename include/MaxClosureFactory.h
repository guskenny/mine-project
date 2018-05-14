#ifndef _MC_FACTORY
#define _MC_FACTORY
#include <stdlib.h>
#include <map>
#include "MaxClosure_Base.h"
#include "graph.h"

class MaxClosureFactory {
  char option; // n=network flow, b,e,p = max flow variants
public:
  MaxClosureFactory(char opt='b');
  ~MaxClosureFactory();
  void setOption(char opt);
  char getOption() const;
  MaxClosure_Base * operator() (const Graph &graph) const;
  MaxClosure_Base * operator() (const Graph &graph, const std::map<int,int> &fixed) const;
};
#endif


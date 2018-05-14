#include "MaxClosureFactory.h"
#include "MaxClosure_BoostMaxFlow_BK.h"
#include "MaxClosure_BoostMaxFlow_EK.h"
#include "MaxClosure_BoostMaxFlow_PR.h"
//#include "MaxClosure_NetworkFlow.h"
#include "MaxClosure_PP.h"
#include <ctype.h>
MaxClosureFactory::MaxClosureFactory(char opt) : option(opt) {}
MaxClosureFactory::~MaxClosureFactory(){}
void MaxClosureFactory::setOption(char opt){ option = opt; }
char MaxClosureFactory::getOption() const {return option; }

MaxClosure_Base *MaxClosureFactory::operator()(const Graph &graph) const
{
      switch(option){
  //    case 'n': return new MaxClosure_NetworkFlow(graph);
      case 'b': return new MaxClosure_BoostMaxFlow_BK(graph);
      case 'e': return new MaxClosure_BoostMaxFlow_EK(graph);
      case 'p': return new MaxClosure_BoostMaxFlow_PR(graph);
      case 'N': case 'B': case 'E': case 'P':
			{ 
				std::map<int,int> fixed;
				MaxClosureFactory mc;
				mc.setOption(tolower(option));
				return new MaxClosure_PP(graph, fixed, mc);
			}
      default:
	  std::cerr << "Unknown MaxClosure option " << option << std::endl;
    }
    return 0;

}

MaxClosure_Base *MaxClosureFactory::operator()(const Graph &graph, const std::map<int,int> &fixed) const
{
      switch(option){
      case 'N': case 'B': case 'E': case 'P': 
      case 'n': case 'b': case 'e': case 'p':
			{ 
				MaxClosureFactory mc;
				mc.setOption(tolower(option));
				return new MaxClosure_PP(graph, fixed, mc);
			}
      default:
	  std::cerr << "Unknown MaxClosure option " << option << std::endl;
    }
    return 0;

}

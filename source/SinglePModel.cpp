#include "SinglePModel.h"
#include <fstream>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>

SinglePModel::SinglePModel(const char *filename) : Daten(filename) {
  defineGraph();
  setUpResourceConsts();
}

SinglePModel::SinglePModel(const Daten &d) : Daten(d) {
  defineGraph();
  setUpResourceConsts();
}

SinglePModel::SinglePModel(const Daten &d, Graph &upitGraph) : Daten(d) {
  defineGraph(upitGraph);
  setUpResourceConsts();
}

// constructor for creating a new SinglePModel after merging
SinglePModel::SinglePModel(SinglePModel *base_model, const std::vector<std::vector<int> > &group_list, std::vector<int> &group_map) : group_map(group_map){

  std::cout << "Creating new aggregate model.. " << std::flush;

  // set daten data
  _rLimits = base_model->getResLimits();
  _t_max = base_model->getNPeriod();
  _nBlock = base_model->getNBlock();
  _r_max = base_model->getnResources();
  _rate = base_model->getDiscountRate();

  // initialise graph to size of group list
  graph = Graph(group_list.size());


  // reset new block group list
  _blocks = std::vector<Block>(group_list.size());

  // create adjacency matrix for new graph so that only one arc gets created
  // between each group pair
  std::vector<std::vector<int> > group_matrix(group_list.size(), std::vector<int> (group_list.size(),0));

  for (int group = 0; group < group_list.size(); ++group){
    // group attributes
    double group_profit = 0;
    std::vector<std::vector<double> > group_rCoef(_t_max, std::vector<double>(_r_max, 0));

    for (int member = 0; member < group_list[group].size(); ++member){
      // get current block id and block object
      int curr_member = group_list[group][member];
      int curr_block_id = curr_member % _nBlock; // gets block id
      int curr_period = curr_member/_nBlock; // gets current period
      const Block & curr_block=base_model->getBlock(curr_block_id);

      // add member to members list of block group
      _blocks[group].addMember(curr_member);

      // add predecessor arcs
      const std::vector<int> pred = curr_block.getPreds();
      int num_pred = curr_block.getNumPred();
      // std::cout << "computing attributes of member " << member << ", block: " << member%_nBlock << ", period: " << member/_nBlock << ", " << num_pred << " preds, pred.size(): " << pred.size() << " " << std::flush;
      for (int p = 0; p < num_pred; ++p){
        // index of previous block
        int prev_block = pred[p];
        // std::cout << prev_block << "/" << std::flush;
        // index of block/time pair
        int prev_member = curr_period * _nBlock + prev_block;
        // std::cout << prev_member%_nBlock << " " << std::flush;
        // add arc to graph if not already there and not the same group
        if (group != group_map[prev_member]){
          if (!group_matrix[group][group_map[prev_member]]){
            // graph.addArc(graph.getNumArcs(), group, group_map[prev_member]);
            _blocks[group].addPred(group_map[prev_member]);
          }
          else{
            group_matrix[group][group_map[prev_member]] = 1;
          }
        }
      }
      // std::cout << "done!" << std::endl;

      // add arc to same block in previous period
      if (curr_period > 0){
        int prev_member = curr_member - _nBlock;
        if (group != group_map[prev_member]){
          if (!group_matrix[group][group_map[prev_member]]){
            // graph.addArc(graph.getNumArcs(), group, group_map[prev_member]);
            _blocks[group].addPred(group_map[prev_member]);
          }
          else{
            group_matrix[group][group_map[prev_member]] = 1;
          }
        }
      }


      // calculate profit for member
      // profit for x[b][t] -> profit for block minus profit for block in next period
      double block_profit=curr_block.getProfit(0);
      if (curr_period < _t_max-1){
        group_profit += block_profit/pow(1+_rate,curr_period) - block_profit/pow(1+_rate,curr_period+1);
      }
      else{
        group_profit += block_profit/pow(1+_rate,_t_max-1);
      }

      // calculate resource usage coefficient for member
      // if a block occurs in the group then the set of time points for that block
      // within the group will be the interval [s, s+t, ..., t]
      // if group is selected, then block is mined at time s, so the resource coefficients
      // for that block will be added to the vector corresponding to period s
      // in order to cancel out subsequent groups, the resource coefficients for the block
      // must be subtracted from the vector corresponding to period t+1, this is because
      // in order for the solution to be feasible, the next group must start at t+1

      // check if current member is the first time the block appears in the group
      if (curr_period == 0 || group != group_map[curr_member - _nBlock]){
        for (int r = 0; r < _r_max; ++r){
          group_rCoef[curr_period][r] += curr_block.getRCoef(0, r);
        }
      }

      // check if current member is the last time the block appears in the group
      if (curr_period == _t_max-1 || group != group_map[curr_member + _nBlock]){
        for (int r = 0; r < _r_max; ++r){
          group_rCoef[curr_period][r] -= curr_block.getRCoef(0, r);
        }
      }
    }

    // set block group profit
    _blocks[group].setProfit(0,group_profit);

    // set block group resource coefficients
    _blocks[group].initialize_rCoef(_t_max, _r_max);
    for (int t = 0; t < _t_max; ++t){
      for (int r = 0; r < _r_max; ++r){
        _blocks[group].setRCoef(t,r,group_rCoef[t][r]);
      }
    }
  }

  std::cout << "done!" << std::endl;

}

SinglePModel::~SinglePModel(){}

int SinglePModel::getnConstraints() const {
  return (int)res.size();
}

double SinglePModel::getProfit(Vertex v) const {
  return profit[v];
}

void SinglePModel::convert_index_to_bd(int index, int &b, int &d){
	const int b_max = getNBlock();
  d = (int) index / b_max;
  b = (int) index % b_max;
}

SinglePModel::Vertex SinglePModel::convert_bd_to_index(int b, int d){
	const int b_max = getNBlock();
	int index = b + b_max * d;
	return index;
}

int SinglePModel::getBlockIdx(SinglePModel::Vertex v) {
  int b,d;
  convert_index_to_bd(v, b, d);
  return b;
}

int SinglePModel::getDestIdx(SinglePModel::Vertex v) {
  int b,d;
  convert_index_to_bd(v, b, d);
  return d;
}

// return vertex at destination 0 for block
SinglePModel::Vertex SinglePModel::lastVertex(int block) {
  return convert_bd_to_index(block,0);
}

// return vertex at destination 0 for block
SinglePModel::Vertex SinglePModel::firstVertex(int block) {
  return convert_bd_to_index(block, getnDestination()-1);
}

// previous destination vertex for same block or InvalidVertex
SinglePModel::Vertex SinglePModel::getPrevDVertex(Vertex v) {
 	const int b = getBlockIdx(v);
	if(v == lastVertex(b) ) return InvalidVertex;
	const int d=getDestIdx(v);
	if( d==0) return InvalidVertex;
	return convert_bd_to_index(b,d-1);
}

// next destination vertex for same block or Invalid vertex
SinglePModel::Vertex SinglePModel::getNextDVertex(Vertex v) {
 	const int b = getBlockIdx(v);
	if(v == firstVertex(b) ) return InvalidVertex;
	const int d=getDestIdx(v);
	if( d+1==getnDestination() ) return InvalidVertex;
	return convert_bd_to_index(b,d+1);
}

void SinglePModel::defineGraph(){
  size_t b_max = getNBlock();
  size_t d_max = getnDestination();
  size_t nv = b_max*d_max;

  profit.resize(nv);

  graph.resizeNodes(nv);

  for (int b=0; b<b_max; ++b){
    const Block &block=getBlock(b);
    // for each block get predecessors and add arcs
	  for(int p=0; p<block.getNumPred(); ++p){
	    const int pred = block.getPreds()[p];
      graph.addArc(pred,b);
	  }
    for (int d=0; d<d_max; ++d){
      if(d>0){
        graph.addArc(convert_bd_to_index(b,d-1), convert_bd_to_index(b,d));
      }
      profit[convert_bd_to_index(b,d)] = block.getProfit(d);
    }
  }
}

void SinglePModel::defineGraph(Graph &upitGraph){
  size_t b_max = getNBlock();
  size_t d_max = getnDestination();
  size_t nv = b_max*d_max;

  profit.resize(nv);

  graph = upitGraph;

  graph.resizeNodes(nv);

  for (int b=0; b<b_max; ++b){
    const Block &block=getBlock(b);
    for (int d=0; d<d_max; ++d){
      if(d>0){
        graph.addArc(convert_bd_to_index(b,d-1), convert_bd_to_index(b,d));
      }
      profit[convert_bd_to_index(b,d)] = block.getProfit(d);
    }
  }
}

// set up resource constraints as a function of vertices
void SinglePModel::setUpResourceConsts(){

  //std::cout << "\ngetnDestination: " << getnDestination() << std::endl
  //          << "getnResources: " << getnResources() << std::endl;

  res.reserve(getnResources());
  const double eps=1e-5;
  for(int r=0;r<getnResources();++r){
    double mult=1;
    if(getResConstrType(r,SINGLE_PERIOD) == 'L'){
	    mult=1;
    }
    else if(getResConstrType(r,SINGLE_PERIOD) == 'G'){
	    mult=-1;
    }
	  else {
	    std::cerr << "WARNING: resource type " << getResConstrType(r,SINGLE_PERIOD)
		  << " not implemented - ignoring constraint " << r << "," << SINGLE_PERIOD
		  << std::endl;
		  continue;
    }
    resLim.push_back(mult*getLimit(r,SINGLE_PERIOD));
    res.resize(res.size()+1);
    auto &resUse = res.back(); // vector of vertex, Rcoef pairs
    for(int b=0;b<getNBlock();++b){
	    int d=getnDestination()-1;		   // start with vertex of last destination
	    Vertex v = convert_bd_to_index(b,d);
	    const Block &block = getBlock(b);
	    // want sum_d Rcoef(d,r)*y_bdt
	    // = sum_d Rcoef(d,r) * (z_bdt - z_earlier(b,d,t))
	    // = sum_d (Rcoef(d,r)-Rcoef(later(b,d,t))) *z_bdt
	    //
//	    Vertex vn = getNextDVertex(v); // get next destination vertex
//	    double q = block.getRCoef(d,r);
//	    if( vn != InvalidVertex && fabs(q) > eps){
//        resUse.push_back(std::make_pair(vn,-mult*q));
//      }
	    while(v!= InvalidVertex){
		    double q = block.getRCoef(getDestIdx(v),r);
		    Vertex vp = getPrevDVertex(v); // get prev destination vertex
		    if( vp != InvalidVertex){
		      q -= block.getRCoef(getDestIdx(vp),r);
        }
		    if( fabs(q) > eps){
		      resUse.push_back(std::make_pair(v,mult*q));
        }
		    v = vp;
	    }
    }
  } // end loop over constraints r
}

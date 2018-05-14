
#include "CumulativeModel.h"
#include <fstream>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
CumulativeModel::CumulativeModel(const char *filename) : Daten(filename) {
	 defineGraph(); 
}
CumulativeModel::CumulativeModel(const Daten &d) : Daten(d) {
	 defineGraph(); 
}
CumulativeModel::~CumulativeModel(){}
int CumulativeModel::getnConstraints() const {
	return (int)res.size();
}
double CumulativeModel::getProfit(int block,int dest,int time) const {
    const Block &b=getBlock(block);
    const double p=b.getProfit(dest);
    switch(time){
	  case 0: return p;
	  case 1: return p/(1+getDiscountRate());
	  case 2: return p/((1+getDiscountRate())*(1+getDiscountRate()));
	  default:
	      return p/pow(1.0+getDiscountRate(),time);
	}
	return 0;
}
double CumulativeModel::getProfit(Vertex v) const {
	return profit[v]; 
}
double CumulativeModel::blockMined(const std::vector<double> &soln,int b,int p) const {
    const double frac = soln[getVertex(b,getnDestination()-1,p)];
    if( p == 0) return frac;
    return frac - soln[getVertex(b,getnDestination()-1,p-1)];
}
double CumulativeModel::blockProcessed(const std::vector<double> &soln,int b,int d,int p) const{
	const Vertex v = getVertex(b,d,p);
	if(v == getEarliestVertex(b)) return soln[v];
	return soln[v] - soln[getEarlierVertex(v)];
}

void CumulativeModel::storeSoln(Sol_Real &sol,const std::vector<double> &y) const
{
    sol.init(getNBlock(),getNPeriod(),getnDestination());
    sol.obj = 0;
    sol.nT = 0;
    for(int b=0;b<getNBlock();++b)
	for(int t=0;t<getNPeriod();++t){
	    sol.x[b][t]=0;
	    for(int d=0;d<getnDestination();++d){
		sol.x[b][t]+= sol.y[b][t][d]
			    = blockProcessed(y,b,d,t);
		sol.obj += sol.y[b][t][d] * getProfit(b,d,t);
	    }
	    if(sol.x[b][t] > 0.01 && t >= sol.nT) sol.nT = t+1;
	}
}
void CumulativeModel::storeSoln(Sol_Int &sol,const std::vector<double> &y) const
{
    sol.init(getNBlock(),getnDestination());
    sol.nT = 0;
    sol.obj=0;
    for(int b=0;b<getNBlock();++b){
	sol.x[b]=getNPeriod();	// never excavated
	for(int t=0;t<getNPeriod();++t){
	    if(blockMined(y,b,t)>0.5){
		sol.x[b]=t;
		if(t>=sol.nT) sol.nT=t+1;
		for(int d=0;d<getnDestination();++d){
		    sol.y[b][d]=blockProcessed(y,b,d,t);
		    sol.obj+=sol.y[b][d]*getProfit(b,d,t);
		}
		break;
	    }
	}
    }
    // normally this shouldn't be necessary, but just in case...
    for(int b=0;b<getNBlock();++b)
	if(sol.x[b] > sol.nT) sol.x[b] = sol.nT; 
    
}
void CumulativeModel::storeSoln(Sol_Int &sol,const std::vector<int> &x) const
{ // convert x to double vector - crude but effective
    std::vector<double> y(x.size());
    for(size_t i=0;i<x.size();++i) y[i] = x[i];
    storeSoln(sol,y);	
}
void CumulativeModel::getSoln(const Sol_Int &sol,std::vector<double> &y) const
{
    const int nBlocks =getNBlock();
    const int d_max =getnDestination();
    y.resize(graph.getNumNodes(),0.0);
    for(auto yy=y.begin();yy!=y.end();++yy) *yy = 0.0; // clear array
    for(int b=0;b<nBlocks;++b)
	if(sol.x[b] < sol.nT){ // block is mined
	    const int t=sol.x[b];
	    double sum=0;
	    for(int d=0;d<d_max;++d){
		sum += sol.y[b][d];
		y[getVertex(b,d,t)] = sum;
	    }
	    for(int v=getVertex(b,t);v!= InvalidVertex;v=getPredVertex(v))
		y[v] = 1;	// block stays mined
	}
}
void CumulativeModel::getSoln(const Sol_Int &sol,std::vector<int> &y) const
{
    std::cerr << "getSoln() not implemented\n";    
}


int CumulativeModel::getBlockIdx(Vertex v) const {
	int b,d,t; 
	convert_index_to_triplet(*this,v,b,t,d); 
	return b; 
}
int CumulativeModel::getDestIdx(Vertex v) const{
	int b,d,t; 
	convert_index_to_triplet(*this,v,b,t,d); 
	return d; 
}
int CumulativeModel::getTimeIdx(Vertex v) const	{
	int b,d,t; 
	convert_index_to_triplet(*this,v,b,t,d); 
	return t; 
}
CumulativeModel::Vertex CumulativeModel::getVertex(int block,int dest,int time) const{ 
	return convert_triplet_to_index(*this,block,time,dest); 
}
CumulativeModel::Vertex CumulativeModel::getSuccVertex(Vertex v) const // successor vertex for same block or InvalidVertex
{ 	const int b = getBlockIdx(v);
	if(v == lastVertex(b) ) return InvalidVertex;
	const int d=getDestIdx(v),t=getTimeIdx(v);
	if( d==0) return getVertex(b,getnDestination()-1,t-1);
	return getVertex(b,d-1,t);
}
CumulativeModel::Vertex CumulativeModel::getPredVertex(Vertex v) const // predecessor vertex (later time) for same block
{ 	const int b = getBlockIdx(v);
	if(v == firstVertex(b) ) return InvalidVertex;
	const int d=getDestIdx(v),t=getTimeIdx(v);
	if( d+1==getnDestination() ) return getVertex(b,0,t+1);
	return getVertex(b,d+1,t);
}
	

double CumulativeModel::calcRedCost(const std::vector<double> &lambda,
				    std::vector<double> &redCost) const
{
    redCost = profit;		// copy basic profit
    double lagConst=0;
    for(size_t r=0;r<resLim.size();++r){
        lagConst += resLim[r]*lambda[r];
	Vertex v; double q;
	BOOST_FOREACH(boost::tie(v,q),res[r]){
	    redCost[v] -= q*lambda[r];
	}
    }
    return lagConst;
}

std::string varName(const CumulativeModel &m,CumulativeModel::Vertex v) {
    return boost::str(boost::format("x%d_b%dd%dt%d")%v%m.getBlockIdx(v)
		      %m.getDestIdx(v)%m.getTimeIdx(v));
}

void CumulativeModel::dump(std::string filename) const
{
	std::ofstream out(filename);
    out << "Problem " << getName() << std::endl;
    out << "Max " << std::endl;
    for(int v=0;v<graph.getNumNodes();++v){
	out << (profit[v] > 0 ? " + " : " - ");
	out << fabs(profit[v]) << " " << varName(*this,v);
	if(v%5==4) out << std::endl;
    }
    out << "Subject To:\n";
    for(long a=0;a<graph.getNumArcs();++a){
	out << (boost::format("A%d: %s - %s >= 0\n")%a
		%varName(*this,graph.getArc(a)->getSrcID())
		%varName(*this,graph.getArc(a)->getTgtID()));
    }
    for(size_t r=0;r<resLim.size();++r){
	out << boost::format("R%d:")%r;
	Vertex v; double q;
	int cnt=0;
	BOOST_FOREACH(boost::tie(v,q),res[r]){
	    if(cnt != 0) out << (q > 0 ? " + " : " - ");
	    out << fabs(q) <<" " << varName(*this,v);
	    if(++cnt % 5==0) out << std::endl;
	}
	out << " <= " << resLim[r] << std::endl;
    }
    out << "Bounds\n";
    for(int v=0;v<graph.getNumNodes();++v){
	if( getPredVertex(v) == InvalidVertex)
	    out << varName(*this,v) << " <= 1\n";
	else if(getSuccVertex(v) == InvalidVertex)
	    out << "0 <= " << varName(*this,v) << std::endl;
    }
    out << "END\n";
}


void CumulativeModel::defineGraph()
{
    size_t nv = getNBlock()*getnDestination()*getNPeriod();
    /*
      blockIdx.resize(nv);
      destIdx.resize(nv);
      timeIdx.resize(nv);
      vertexIdx.resize(getNBlock());
    int v=0;
    for(int b=0;b<getNBlock();++b){
	vertexIdx[b].resize(getnDestination());
	for(int d=0;d<getnDestination();++d){
	    vertexIdx[b][d].resize(getNPeriod());
	    for(int t=0;t<getNPeriod();++t){
		vertexIdx[b][d][t]=v;
		blockIdx[v]=b;
		destIdx[v] =d;
		timeIdx[v] =t;
		++v;
	    }
	}
    } 
    */
    profit.resize(nv);
    // define profit for each vertex
    for(Vertex v=0;v<nv;++v){
	const int b=getBlockIdx(v),d=getDestIdx(v),t=getTimeIdx(v);
	profit[v] = getProfit(b,d,t);
	Vertex vv=getLaterVertex(v);
	if(vv != InvalidVertex)
	    profit[v] -= getProfit(getBlockIdx(vv),getDestIdx(vv),getTimeIdx(vv));
    }

    // define the actual graph
    graph.resizeNodes(nv);
    for(int b=0;b<getNBlock();++b){
	const Block &block=getBlock(b);
	for(int p=0;p<block.getNumPred();++p){
	    const int pred = block.getPreds()[p];
	    const int d=getnDestination()-1;
	    for(int t=0;t<getNPeriod();++t){
		graph.addArc(getVertex(pred,d,t),getVertex(b,d,t));
	    }
	}
	for(Vertex v=lastVertex(b);v!=firstVertex(b);v=getPredVertex(v))
	    graph.addArc(getPredVertex(v),v);
    }
  
    // set up resource constraints as a function of vertices
    res.reserve(getnResources()*getNPeriod());
    const double eps=1e-5;
    for(int t=0;t<getNPeriod();++t)
      for(int r=0;r<getnResources();++r){
	    double mult=1;
	    if(getResConstrType(r,t) == 'L')
			mult=1;
	    else if(getResConstrType(r,t) == 'G')
			mult=-1;
			else{
			  std::cerr << "WARNING: resource type " << getResConstrType(r,t)
				 << " not implemented - ignoring constraint " << r << "," << t
				 << std::endl;
				continue;
		}
	    resLim.push_back(mult*getLimit(r,t));
	    res.resize(res.size()+1);
	    auto &resUse = res.back(); // vector of vertex, Rcoef pairs
	    for(int b=0;b<getNBlock();++b){
			int d=0;		   // start with earliest vertex
			Vertex v=getVertex(b,d,t);	
			const Block &block=getBlock(b);
			// want sum_d Rcoef(d,r)*y_bdt
			// = sum_d Rcoef(d,r) * (z_bdt - z_earlier(b,d,t))
			// = sum_d (Rcoef(d,r)-Rcoef(later(b,d,t))) *z_bdt
			//
			Vertex ve = getEarlierVertex(v);
			double q = block.getRCoef(d,r);
			if( ve != InvalidVertex && fabs(q) > eps) resUse.push_back(std::make_pair(ve,-mult*q));
			while(v!= InvalidVertex && getTimeIdx(v) == t){
				q = block.getRCoef(getDestIdx(v),r);
				Vertex vl = getLaterVertex(v);
				if( vl != InvalidVertex && getTimeIdx(vl) == t)
				q -= block.getRCoef(getDestIdx(vl),r);
				if( fabs(q) > eps)
				resUse.push_back(std::make_pair(v,mult*q));
				v = vl;					 
			}
		}
	} // end loop over constraints r,t
} // end defineGraph()


// extract information from BranchNode info struct
void CumulativeModel::BN_info_to_map(const BranchNode_info &info, std::map<int,int> &fixed){
    fixed.clear();
  for (size_t b = 0; b < info.time.size(); ++b){
    if (info.time[b][AFTER] - info.time[b][BEFORE] <= 0){ // never mine
	fixed[getLatestVertex(b)] = 0;	// still not mined at latest time
    } else { // must be mined at specific time or within interval
      if (info.time[b][BEFORE] > 0){
	  fixed[getVertex(b, info.time[b][BEFORE] - 1)] = 0;
      }
      if( info.time[b][AFTER] < getNPeriod() )
	  fixed[getVertex(b, info.time[b][AFTER]-1)] = 1;
    }
  }  
}



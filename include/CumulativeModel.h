/*
 * This class defines the precedence constrained pit scheduling problem
 * (PCPSP) in terms of a cumulative model (that is with variables representing
 * the total amount mined for a block up to that time/destination)
 * The fromulation defiens a graph of precedences with nodes representing
 * (block,destination,time) triplets, profit per node and resource
 * constraints (<= only!)
 * Note that the decision variables associated with a node in this formulation
 * correspond to the cumulative sum of the total fraction of a block b that
 * has been sent to destination d (or earlier) by time t (or earlier)
 */
#ifndef __CumulativeModel_H__
#define __CumulativeModel_H__
#include "daten.h"
#include "graph.h"
#include <map>
#include "util.h"
#include <math.h>
#include "BranchNode.h"

#define InvalidVertex -1
#define BEFORE 0
#define AFTER 1

/// CumulativeModel : defined in terms of variables z_bdt such that
/// z_bdt = fraction of block b mined and processed up to d , t
///       = sum_t'<t  sum_d' y_bd't' + sum_d'<=d y_bd't
/// ie ordering by increasing time and destination index
/// Hence y_bdt = z_bdt - z_earlier(b,d,t)
class CumulativeModel : public Daten {
public:
    typedef int Vertex;
    //const int InvalidVertex;
    // constructors - load from file or copy from previously loaded
    CumulativeModel(const char *filename); // load data or
    CumulativeModel(const Daten &d);
	~CumulativeModel();
	
    int getnConstraints() const;
    double getProfit(int block,int dest,int time) const;
    // profit of a vertex = profit(v)-profit(later(v))
    double getProfit(Vertex v) const;
    const std::vector<double> &getProfit() const { return profit; }
    // calculate reduced cost and return lagrangian constant
    // assumes lambda >= 0 is length getnConstraints()
    double calcRedCost(const std::vector<double> &lambda,std::vector<double> &redCost) const;
    // more or less the same as evaluate_supply(data,lambda,supply) in util.h
    
    // inspecting a solution:
    // blockMined(s,b,t) = fraction of block b mined during period t
    double blockMined(const std::vector<double> &soln,int b,int p) const;
    // blockProcessed(s,b,d,t) fraction of b sent to destination d in period t
    double blockProcessed(const std::vector<double> &soln,int b,int d,int p) const;
    // the following methods are really static
    void storeSoln(Sol_Real &sol,const std::vector<double> &y) const;
    void storeSoln(Sol_Int &sol,const std::vector<double> &y) const;
    void storeSoln(Sol_Int &sol,const std::vector<int> &y) const;
    void getSoln(const Sol_Int &sol,std::vector<double> &y) const;
    void getSoln(const Sol_Int &sol,std::vector<int> &y) const;

    
    // indexing between block,destination,period triples and vertices    
    int getBlockIdx(Vertex v) const;
    int getDestIdx(Vertex v) const;
    int getTimeIdx(Vertex v) const;
    Vertex getVertex(int block, int dest,int time) const;
    Vertex getVertex(int block,int time) const
        // vertex corresponding to the  binary variable block,time variable
	{ return getVertex(block,getnDestination()-1,time); }

    // extract fixed information from BranchNode_info struct
    void BN_info_to_map(const BranchNode_info &info, std::map<int,int> &fixed);

    //------------------------------------------------------------------
    // Vertices for the same block form chains of successive vertices
    // Note that successor represents earlier in time.
    Vertex getSuccVertex(Vertex v) const; // successor vertex for same block or InvalidVertex
    Vertex getPredVertex(Vertex v) const; // predecessor vertex (later time) for same block
    Vertex getEarlierVertex(Vertex v) const {return getSuccVertex(v); }
    Vertex getLaterVertex(Vertex v) const {return getPredVertex(v);   }
    Vertex getEarliestVertex(int block) const {return lastVertex(block);  }
    Vertex getLatestVertex(int block) const {return firstVertex(block);  }
    Vertex firstVertex(int block) const // vertex without predecessor (end time)
	{ return getVertex(block,getnDestination()-1,getNPeriod()-1); }
    Vertex lastVertex(int block) const // vertex without successort (init time)
	{ return getVertex(block,0,0); }
    void dump(std::string filename) const; // write problem out to file in LP format


 
// Problem information - should be treated as read-only
    Graph graph;
    // resource constraints for x[v] = fraction of vertex v
    //  sum_i  res[r][i].second * x[res[r][i].first] <= resLim[r]
    std::vector<double> resLim; // resource limit (RHS)
    std::vector<std::vector<std::pair<Vertex,double> > > res; 
  
  
protected:
    std::vector<double> profit; // profit for each vertex    
    void defineGraph(); // create cumulative graph and resource constraints
    // std::vector<int> blockIdx;
    // std::vector<int> destIdx;
    // std::vector<int> timeIdx;
    // std::vector<std::vector<std::vector<Vertex> > > vertexIdx;

};
    
#endif

#ifndef __SinglePModel_H__
#define __SinglePModel_H__
#include "daten.h"
#include "graph.h"
#include "util.h"
#include "UpitSolver.h"
#include <math.h>

#define InvalidVertex -1
#define SINGLE_PERIOD 0


class SinglePModel : public Daten {
  public:
    typedef int Vertex;
    SinglePModel(const char *filename);
    SinglePModel(const Daten &d);
    SinglePModel(const Daten &d, Graph &upitGraph);
    SinglePModel(SinglePModel *base_model, const std::vector<std::vector<int> > &group_list, std::vector<int> &group_map);
    ~SinglePModel();

    int getnConstraints() const;
    const std::vector<double> &getProfit() const { return profit; };
    double getProfit(Vertex v) const;
    void convert_index_to_bd(int index, int &b, int &d);
    Vertex convert_bd_to_index(int b, int d);
    int getBlockIdx(Vertex v);
    int getDestIdx(Vertex v);
    Vertex lastVertex(int block);
    Vertex firstVertex(int block);
    Vertex getPrevDVertex(Vertex v);
    Vertex getNextDVertex(Vertex v);

    Graph graph;

    std::vector<double> resLim; // resource limit (RHS)
    std::vector<std::vector<std::pair<Vertex,double> > > res;
    std::vector<int> group_map; // used for merged group model

  protected:
    std::vector<double> profit;
    void defineGraph();
    void defineGraph(Graph &upitGraph);
    void setUpResourceConsts();
};

#endif

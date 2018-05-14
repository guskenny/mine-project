TARGET = lraco_opbs
INC = /usr/local/gurobi/6.5.1/include/
LIB = /usr/local/gurobi/6.5.1/lib/ 
CPLEX_ROOT = /usr/local/cplex/12.6.3/
LIB_CPLEX = $(CPLEX_ROOT)cplex/lib/x86-64_linux/static_pic/
LIB_CONCERT = $(CPLEX_ROOT)concert/lib/x86-64_linux/static_pic/
CXX = mpic++
CXXFLAGS = -m64 -O3 -DIL_STD -fopenmp -std=c++0x -pthread -wd873 -g 
CPPFLAGS = -I ../include -I$(INC) -I$(CPLEX_ROOT)concert/include/ -I$(CPLEX_ROOT)cplex/include/ -std=c++0x -pthread
LDFLAGS = -L${LIB} -L${LIB_CPLEX} -L${LIB_CONCERT}  
LIBS = -L${LIB}$ -L${LIB_CPLEX} -L${LIB_CONCERT} -lilocplex -lcplex -lconcert -lgurobi_c++ -lpthread -lm -lgurobi65
OBJS = Random.o Timer.o ../source/daten.o pheromones.o mip.o ../source/lr_graph.o ../source/graph.o ../source/MaxClosure_Base.cpp MaxClosure_NetworkFlow_LR.o solver.o solver_functions.o ACO_solution.o

all: ${TARGET}

lraco_opbs : lraco_opbs.cpp $(OBJS)
	${CXX} ${CXXFLAGS} ${CPPFLAGS} $^ -o $@  ${LIBS}
clean:
	rm -f *.o ../source/*.o *~ ${TARGET} core gmon.out


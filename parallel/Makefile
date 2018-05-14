All:
	cd ../;make -j30

CPPFLAGS = -std=c++0x -pthread -wd873
# Separate compile options per configuration
ifeq ($(CONFIG),debug)
CPPFLAGS += -g -O0
else
CPPFLAGS += -O3
endif

CPP=mpic++
CC=mpicc

INC = /usr/local/gurobi/6.5.1/include/
LIB = /usr/local/gurobi/6.5.1/lib/
CPLEX_ROOT = /usr/local/cplex/12.6.3/
LIB_CPLEX = $(CPLEX_ROOT)cplex/lib/x86-64_linux/static_pic/
LIB_CONCERT = $(CPLEX_ROOT)concert/lib/x86-64_linux/static_pic/

GCC_VER=$(shell $(CXX) -dumpversion|grep -o '[0-9].[0-9][0-9]*')
CXXFLAGS = -m64 -DIL_STD -fopenmp
CPPFLAGS += -I../include -I$(INC) -I$(CPLEX_ROOT)concert/include/ -I$(CPLEX_ROOT)cplex/include/
LDFLAGS = -L${LIB} -L${LIB_CPLEX} -L${LIB_CONCERT}  
LIBS = -L${LIB}$ -L${LIB_CPLEX} -L${LIB_CONCERT} -lilocplex -lcplex -lconcert -lgurobi_c++ -lpthread -lm -lgurobi65

LAPSO_SRC=../source/Particle.cpp $(wildcard QOL/lagrangian/[^d]*.cpp)
LR = ../LagrangianHeuristic

SRC := main.cpp ../source/daten.cpp ../source/network.cpp ../source/boostMaxFlow.cpp ../source/CumulativeModel.cpp ../source/MaxClosure_NetworkFlow.cpp ../source/util.cpp ../source/MaxClosure_PP.cpp $(wildcard ../source/MaxClosure_B*.cpp) ../source/MaxClosureFactory.cpp ../source/UpitSolver.cpp ../source/graph.cpp $(LR)/solver.cpp $(LR)/mip.cpp $(LR)/ACO_solution.cpp $(LR)/Timer.cpp $(LR)/Random.cpp $(LR)/MaxClosure_NetworkFlow_LR.cpp $(LR)/pheromones.cpp ../source/lr_graph.cpp $(LR)/solver_functions.cpp ../source/solver_interface.cpp ../source/bz.cpp ../source/rlp.cpp ../source/BranchBound.cpp MineProblem.cpp $(LAPSO_SRC)

INCS := $(wildcard *.h)
#INC := $(SRC:%.cpp=%.h)
OBJ := $(SRC:%.cpp=%.o)

PROGRAM = test

default: $(PROGRAM)

#Rebuild *.cpp
$(OBJ): %.o : %.cpp $(INCS)
	$(CPP) ${CXXFLAGS} ${CPPFLAGS} -c $< -o $@

$(PROGRAM): $(OBJ)
	$(CPP) ${CXXFLAGS} ${CPPFLAGS} -o $(PROGRAM) $(OBJ) $(LIBS) $(LDFLAGS)

clean:
	/bin/rm -f *.o $(PROGRAM)


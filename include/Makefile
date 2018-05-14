SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEX_ROOT = /usr/local/cplex/12.6.3/
LIB_CPLEX = $(CPLEX_ROOT)cplex/lib/x86-64_linux/static_pic/
LIB_CONCERT = $(CPLEX_ROOT)concert/lib/x86-64_linux/static_pic/
#CXXFLAGS = -m64 -O3 -DIL_STD -fopenmp
#CPPFLAGS = -I ../include -I$(INC) -I$(CPLEX_ROOT)concert/include/ -I$(CPLEX_ROOT)cplex/include/

CPLEXDIR      = /usr/local/cplex/12.6.3/cplex/lib/x86-64_linux/static_pic/
CONCERTDIR    = /usr/local/cplex/12.6.3/concert/lib/x86-64_linux/static_pic/
# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

#CCC = g++ -O0
#CC  = gcc -O0
CCC = g++ 
CC  = gcc

# user working directory
USERDIR       = /home/dthiruva/BHPB-MAXIMA/Code/OPBS/
USERSCRDIR    = $(USERDIR)/source
USERINCDIR    = $(USERDIR)/include
USERDATA      = $(USERDIR)/data

# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------
#CCOPT = -m64 -g -fPIC -fno-strict-aliasing -fexceptions -DIL_STD
CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD
COPT  = -m64 -fPIC -fno-strict-aliasing


# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) -L${LIB_CPLEX} -L${LIB_CONCERT}  
CLNDIRS   = -L$(CPLEXLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread
CLNFLAGS  = -lcplex -lm -lpthread


all:
	make all_cpp

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(USERINCDIR) -I$(CPLEX_ROOT)concert/include/ -I$(CPLEX_ROOT)cplex/include/ -I../source/
CFLAGS  = $(COPT)  -I$(CPLEXINCDIR)


#------------------------------------------------------------
#  make all      : to compile the examples. 
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------
CPP_EX = LPR_pcpsp

all_cpp: $(CPP_EX)

# ------------------------------------------------------------

clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm -rf $(C_EX) $(CPP_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp

# ------------------------------------------------------------
LPR_pcpsp: daten.o util.o max_closure.o rlp.o lpr_bz.o 
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o LPR_pcpsp daten.o util.o max_closure.o rlp.o lpr_bz.o $(CCLNFLAGS)

daten.o: $(USERSCRDIR)/daten.cpp 
	$(CCC) $(CCFLAGS) -c $(USERSCRDIR)/daten.cpp
	
#check.o: $(USERSCRDIR)/check.cpp 
#	$(CCC) $(CCFLAGS) -c $(USERSCRDIR)/check.cpp
	
util.o: $(USERSCRDIR)/util.cpp 
	$(CCC) $(CCFLAGS) -c $(USERSCRDIR)/util.cpp 

max_closure.o: $(USERSCRDIR)/max_closure.cpp 
	$(CCC) $(CCFLAGS) -c $(USERSCRDIR)/max_closure.cpp

rlp.o: $(USERSCRDIR)/rlp.cpp 
	$(CCC) $(CCFLAGS) -c $(USERSCRDIR)/rlp.cpp

lpr_bz.o: $(USERSCRDIR)/lpr_bz.cpp
	$(CCC) $(CCFLAGS) -c $(USERSCRDIR)/lpr_bz.cpp 

	


# Local Variables:
# mode: makefile
# End:

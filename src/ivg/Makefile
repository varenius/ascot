all: ivgUtilsLib

DEFINES=
#DEFINES= -DDEBUG_REFFRAME=2 -DDEBUG_VLBI=2 -DDEBUG_SOLVER=2 -DDEBUG_TROP=1
#DEFINES= -DDEBUG_VLBI=2 -DDEBUG_SOLVER=1 -DDEBUG_TROP=1 
#DEFINES= -DDEBUG_TROP=1

CXX=g++
CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} -DBOOST_MATH_DISABLE_FLOAT128
CXXOPTS=-O3 -std=c++14 $(DEFINES) -g -fPIC

LIBDIR = $(ASCOT_DIR)/lib
IVGUTILS = .


INCBLAS  = -I$(OB_DIR)/include

INCALL = $(INCBLAS) -I./ \
          -I$(ASCOT_DIR)/include \
	  -I$(ASCOT_DIR)/src/projectpluto/include \
	  -I$(ASCOT_DIR)/src/iers\
	  -I$(ASCOT_DIR)/src/sofa/include
	
# lpsked needs qt
INCQT = -I. -I/usr/include/x86_64-linux-gnu/qt5/QtCore -I/usr/include/x86_64-linux-gnu/qt5/QtGui -I/usr/include/x86_64-linux-gnu/qt5/


OBJS = matrix.o funother.o date.o auxfunc.o tictoc.o sparse2array.o eop_series.o station.o analysis_station.o parser.o trf.o session.o scan.o obs.o troposphere.o turbulence.o\
       source.o crf.o param.o param_list.o logger.o ls_solution.o lsa.o lsc.o icls.o session_inout.o ls_neq.o vgosdb.o simulation.o schedule.o sinex.o masterfile.o db_download.o\
       wrapper.o ltnStation.o log2nc.o ltnSession.o cableCor.o atmosphereParam.o transits.o grid.o milp.o hemisphereData.o

matrix.o: $(IVGUTILS)/matrix.cpp $(IVGUTILS)/matrix.h
	g++ $(OPTIMIZE) $(DEFINES) -c $(IVGUTILS)/matrix.cpp $(INCALL)

funother.o: $(IVGUTILS)/funother.cpp $(IVGUTILS)/funother.h
	g++ $(OPTIMIZE) -c $(IVGUTILS)/funother.cpp -I./

date.o: $(IVGUTILS)/date.cpp
	$(CXX) $(CXXOPTS) -c $(IVGUTILS)/date.cpp $(OPTIMIZE) $(DEFINES) $(INCALL)

db_download.o: db_download.cpp db_download.h
	$(CXX) $(CXXOPTS) -c db_download.cpp $(INCALL) -I$(ASCOT_DIR)/src/ltn

eop_series.o: $(IVGUTILS)/eop_series.cpp $(IVGUTILS)/eop_series.h
	$(CXX) $(CXXOPTS) -c $(IVGUTILS)/eop_series.cpp $(OPTIMIZE) $(DEFINES) $(INCALL)

station.o: station.cpp station.h
	$(CXX) $(CXXOPTS) -c station.cpp $(INCALL) 

analysis_station.o: analysis_station.cpp analysis_station.h
	$(CXX) $(CXXOPTS) -c analysis_station.cpp $(INCALL) 

auxfunc.o: $(IVGUTILS)/auxfunc.cpp  $(IVGUTILS)/auxfunc.h
	$(CXX) $(CXXOPTS) -o auxfunc.o -c $(IVGUTILS)/auxfunc.cpp $(INCALL)

session.o: session.cpp session.h
	$(CXX) $(CXXOPTS) -c session.cpp $(INCALL)  $(INCQT)

scan.o: scan.cpp scan.h
	$(CXX) $(CXXOPTS) -c scan.cpp $(INCALL) 

obs.o: obs.cpp obs.h
	$(CXX) $(CXXOPTS) -c obs.cpp $(INCALL)  $(INCQT)

troposphere.o: troposphere.cpp troposphere.h
	$(CXX) $(CXXOPTS) -c troposphere.cpp $(INCALL) 
	
turbulence.o: turbulence.cpp turbulence.h
	$(CXX) $(CXXOPTS) -c turbulence.cpp $(INCALL) #-fopenmp

parser.o: parser.h parser.cpp
	$(CXX) $(CXXOPTS) -c parser.cpp $(INCALL) 
	
trf.o: trf.h trf.cpp
	$(CXX) $(CXXOPTS) -c trf.cpp $(INCALL) 
	
crf.o: crf.h crf.cpp
	$(CXX) $(CXXOPTS) -c crf.cpp $(INCALL) 
	
source.o: source.h source.cpp
	$(CXX) $(CXXOPTS) -c source.cpp $(INCALL) 
		
param.o: param.h param.cpp
	$(CXX) $(CXXOPTS) -c param.cpp $(INCALL) 
	
param_list.o: param_list.h param_list.cpp
	$(CXX) $(CXXOPTS) -c param_list.cpp $(INCALL) 

session_inout.o: session_inout.h session_inout.cpp
	$(CXX) $(CXXOPTS) -c session_inout.cpp $(INCALL)  $(INCQT)

vgosdb.o: vgosdb.h vgosdb.cpp
	$(CXX) $(CXXOPTS) -c vgosdb.cpp $(INCALL)
	
tictoc.o:  $(IVGUTILS)/tictoc.h
	g++ -c -o tictoc.o tictoc.h

logger.o: $(IVGUTILS)/logger.h logger.cpp
	$(CXX) $(CXXOPTS) -c logger.cpp $(INCALL)

lsm.o: lsm.h lsm.cpp
	$(CXX) $(CXXOPTS) -c lsm.cpp $(INCALL) 

ls_solution.o: ls_solution.h ls_solution.cpp
	$(CXX) $(CXXOPTS) -c ls_solution.cpp $(INCALL) 

ls_neq.o: ls_neq.h ls_neq.cpp
	$(CXX) $(CXXOPTS) -c ls_neq.cpp $(INCALL) 

lsa.o: lsa.h lsa.cpp
	$(CXX) $(CXXOPTS) -c lsa.cpp $(INCALL) 

lsc.o: lsc.h lsc.cpp
	$(CXX) $(CXXOPTS) -c lsc.cpp $(INCALL) 
	
icls.o: icls.h icls.cpp
	$(CXX) $(CXXOPTS) -c icls.cpp $(INCALL) 

sparse2array.o: $(IVGUTILS)/sparse2array.c
	   gcc -O3 -c $(IVGUTILS)/sparse2array.c
	
hemisphereData.o: hemisphereData.h hemisphereData.cpp
	$(CXX) $(CXXOPTS) -c hemisphereData.cpp $(INCQT) $(INCALL)	

simulation.o: simulation.h simulation.cpp
	$(CXX) $(CXXOPTS) -c simulation.cpp $(INCALL)  $(INCQT)
	
schedule.o: schedule.h schedule.cpp
	$(CXX) $(CXXOPTS) -c schedule.cpp $(INCALL)  $(INCQT)
	
sinex.o: sinex.h sinex.cpp
	$(CXX) $(CXXOPTS) -c sinex.cpp $(INCALL) 
	
masterfile.o: masterfile.h masterfile.cpp
	$(CXX) $(CXXOPTS) -c masterfile.cpp $(INCALL)
	
wrapper.o: wrapper.h wrapper.cpp
	$(CXX) $(CXXOPTS) -c wrapper.cpp $(INCALL)
	
ltnStation.o: ltn/ltnStation.cpp ltn/ltnStation.h
	$(CXX) $(CXXOPTS)  -c ltn/ltnStation.cpp $(INCALL)
	
log2nc.o: ltn/log2nc.cpp ltn/log2nc.h
	$(CXX) $(CXXOPTS)  -c ltn/log2nc.cpp $(INCALL)
	
ltnSession.o: ltn/ltnSession.cpp ltn/ltnSession.h
	$(CXX) $(CXXOPTS)  -c ltn/ltnSession.cpp  $(INCALL)
		
cableCor.o: ltn/cableCor.cpp ltn/cableCor.h
	$(CXX) $(CXXOPTS)  -c ltn/cableCor.cpp $(INCALL)
	
atmosphereParam.o: ltn/atmosphereParam.cpp ltn/atmosphereParam.h
	$(CXX) $(CXXOPTS)  -c ltn/atmosphereParam.cpp $(INCALL)
	
transits.o: lpSked/transits.cpp lpSked/transits.h
	$(CXX) $(CXXOPTS) -c lpSked/transits.cpp $(INCQT) $(INCALL)
	
grid.o: lpSked/grid.cpp lpSked/grid.h
	$(CXX) $(CXXOPTS) -c lpSked/grid.cpp $(INCALL) $(INCQT)
	
milp.o: lpSked/milp.cpp lpSked/milp.h
	$(CXX) $(CXXOPTS) -c lpSked/milp.cpp $(INCALL) $(INCQT) -I $(GUROBI_HOME)/include

ivgUtilsLib: $(OBJS)
	rm -f $(LIBDIR)/libIvgUtils.a
	ar rs $(LIBDIR)/libIvgUtils.a $(OBJS)

clean: 
	rm -f $(OBJS)
	rm -f $(LIBDIR)/libIvgUtils.a

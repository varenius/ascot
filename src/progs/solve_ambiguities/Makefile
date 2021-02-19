#solve_ambiguities makefile
USE_AHC = FALSE

CXX=g++

ifeq ($(USE_AHC),TRUE)
CXXOPTS=-O3 -std=c++14 -g -DUSE_AHC -DBOOST_MATH_DISABLE_FLOAT128 #-DDEBUG
else
CXXOPTS=-O3 -std=c++14 -g -fPIC -DBOOST_MATH_DISABLE_FLOAT128 #-DDEBUG
endif


LIB_DIR = $(ASCOT_DIR)/lib
BIN_DIR = $(ASCOT_DIR)/bin

INCBLAS  = -I$(OB_DIR)/include
LIBBLAS  = -L$(OB_DIR)/lib/ -llapack -lopenblas  -lpthread -lgfortran

INCQT = -I. -I/usr/include -I/usr/include/x86_64-linux-gnu/qt5/QtGui/ -I/usr/include/x86_64-linux-gnu/qt5/QtCore -I/usr/include/x86_64-linux-gnu/qt5/QtWidgets -I/usr/include/x86_64-linux-gnu/qt5/QtPrintSupport
LIBQT = -lQt5Gui -L/usr/lib64 -L/usr/X11R6/lib -lQt5Core -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lQt5Widgets -lQt5PrintSupport -lboost_program_options


ifeq ($(USE_AHC),TRUE)
ALGLIB_DIR = $(ASCOT_DIR)/src/alglib
INC_ALL = -I./ $(INCQT) \
			 $(INCBLAS) \
			 -I$(ASCOT_DIR)/src/iers \
			 -I$(ASCOT_DIR)/src/ivg \
			 -I$(ASCOT_DIR)/include \
			 -I$(ASCOT_DIR)/src/projectpluto/include \
			 -I$(ASCOT_DIR)/src/sofa/include \
			 -I$(ASCOT_DIR)/src/qt_plot \
			 -I$(ASCOT_DIR)/src/analysistools \
			 -I$(ASCOT_DIR)/src/ltn \
			 -I$(ALGLIB_DIR)   
else
INC_ALL = -I./ $(INCQT) \
			 $(INCBLAS) \
			 -I$(ASCOT_DIR)/src/iers \
			 -I$(ASCOT_DIR)/src/ivg \
			 -I$(ASCOT_DIR)/include \
			 -I$(ASCOT_DIR)/src/projectpluto/include \
			 -I$(ASCOT_DIR)/src/sofa/include \
			 -I$(ASCOT_DIR)/src/qt_plot \
			 -I$(ASCOT_DIR)/src/analysistools \
			 -I$(ASCOT_DIR)/src/ltn
endif


IGG_VLBI_LIB = -L$(ASCOT_DIR)/lib -lIvgUtils -lqt_plot -lAnalysisTools $(ASCOT_DIR)/src/iers/libIers2010.a

LIBS =						\
	-lcurl -lconfig++ \
	$(IGG_VLBI_LIB) \
	-L$(ASCOT_DIR)/src/projectpluto/lib -ljpleph -lsat\
	-L$(ASCOT_DIR)/src/sofa/lib/ -lsofa_c \
	$(LIBBLAS) \
	$(LIBLAPACK) -lmatio -lhdf5 -L/usr/include/mysql -lmysqlclient \
	$(LIBQT) \
	-lgmp \
	-lboost_regex \
	-lnetcdf_c++ -lnetcdf -lfftw3 `gsl-config --libs`


ifeq ($(USE_AHC),TRUE)
OBJS = network.o triangle.o alglibinternal.o alglibmisc.o ap.o dataanalysis.o optimization.o linalg.o statistics.o specialfunctions.o solvers.o
else
OBJS = network.o triangle.o
endif

all: $(OBJS) solve_ambiguities
	
network.o: network.cpp network.h
	$(CXX) $(CXXOPTS)  -c network.cpp $(LIBS) $(INC_ALL)
	
triangle.o: triangle.cpp triangle.h
	$(CXX) $(CXXOPTS)  -c triangle.cpp $(LIBS) $(INC_ALL)

ifeq ($(USE_AHC),TRUE)
alglibinternal.o: $(ALGLIB_DIR)/alglibinternal.cpp $(ALGLIB_DIR)/alglibinternal.h
	$(CXX) $(CXXOPTS)  -c $(ALGLIB_DIR)/alglibinternal.cpp $(LIBS) $(INC_ALL) 
	
alglibmisc.o: $(ALGLIB_DIR)/alglibmisc.cpp $(ALGLIB_DIR)/alglibmisc.h
	$(CXX) $(CXXOPTS)  -c $(ALGLIB_DIR)/alglibmisc.cpp $(LIBS) $(INC_ALL) 
	
ap.o: $(ALGLIB_DIR)/ap.cpp $(ALGLIB_DIR)/ap.h
	$(CXX) $(CXXOPTS)  -c $(ALGLIB_DIR)/ap.cpp $(LIBS) $(INC_ALL) 
	
dataanalysis.o: $(ALGLIB_DIR)/dataanalysis.cpp $(ALGLIB_DIR)/dataanalysis.h
	$(CXX) $(CXXOPTS)  -c $(ALGLIB_DIR)/dataanalysis.cpp $(LIBS) $(INC_ALL) 
	
linalg.o: $(ALGLIB_DIR)/linalg.cpp $(ALGLIB_DIR)/linalg.h
	$(CXX) $(CXXOPTS)  -c $(ALGLIB_DIR)/linalg.cpp $(LIBS) $(INC_ALL)
	
optimization.o: $(ALGLIB_DIR)/optimization.cpp $(ALGLIB_DIR)/optimization.h
	$(CXX) $(CXXOPTS)  -c $(ALGLIB_DIR)/optimization.cpp $(LIBS) $(INC_ALL)
	
solvers.o: $(ALGLIB_DIR)/solvers.cpp $(ALGLIB_DIR)/solvers.h
	$(CXX) $(CXXOPTS)  -c $(ALGLIB_DIR)/solvers.cpp $(LIBS) $(INC_ALL)
	
specialfunctions.o: $(ALGLIB_DIR)/specialfunctions.cpp $(ALGLIB_DIR)/specialfunctions.h
	$(CXX) $(CXXOPTS)  -c $(ALGLIB_DIR)/specialfunctions.cpp $(LIBS) $(INC_ALL)
	
statistics.o: $(ALGLIB_DIR)/statistics.cpp $(ALGLIB_DIR)/statistics.h
	$(CXX) $(CXXOPTS)  -c $(ALGLIB_DIR)/statistics.cpp $(LIBS) $(INC_ALL)
endif
	
solve_ambiguities: solve_ambiguities.cpp
	$(CXX) $(CXXOPTS) -o $(BIN_DIR)/solve_ambiguities solve_ambiguities.cpp  $(OBJS)  $(LIBS) $(INC_ALL) 
	
clean:
	rm -f $(BIN_DIR)/solve_ambiguities $(OBJS)








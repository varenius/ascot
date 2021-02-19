USEMPI=FALSE


ifeq ($(USEMPI),TRUE)
    DEFINES=-DQT_NO_DEBUG -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED
    CXX = $(MPI_DIR)/bin/mpic++
    CXXOPTS = -std=c++14 -O3 -DUSE_MPI -g -DBOOST_MATH_DISABLE_FLOAT128
    INCMPI = -I$(MPI_DIR)/include
else
    DEFINES=-DQT_NO_DEBUG -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED -DUSE_MPI
    CXXOPTS=-O3 -std=c++14 -g -fPIC -DBOOST_MATH_DISABLE_FLOAT128
    CXX=g++
    INCMPI =
endif


LIB_DIR = $(ASCOT_DIR)/lib
BIN_DIR = $(ASCOT_DIR)/bin

INCBLAS  = -I$(OB_DIR)/include
LIBBLAS  = -L$(OB_DIR)/lib/ -llapack -lopenblas -lgfortran -lpthread

INCQT = -I. -I/usr/include/x86_64-linux-gnu/qt5/QtGui -I/usr/include/x86_64-linux-gnu/qt5/QtCore -I/usr/include/x86_64-linux-gnu/qt5/QtWidgets -I/usr/include/x86_64-linux-gnu/qt5/QtPrintSupport -I/opt/gurobi911/linux64/include
LIBQT = -lQt5Gui -L/usr/lib64 -L/usr/X11R6/lib -lQt5Core -L/usr/lib/x86_64-linux-gnu/hdf5/serial -L/opt/gurobi911/linux64/lib -lQt5Widgets -lQt5PrintSupport -lboost_program_options 

#INC_ALL = -I./ \
#			 $(INCBLAS) \
#			 -I$(ASCOT_DIR)/src/ivg \
#			 -I$(ASCOT_DIR)/include \
#			 -I$(ASCOT_DIR)/src/projectpluto/include \
#			 -I$(ASCOT_DIR)/src/sofa/include

INC_ALL = -I./ $(INCQT) \
			 $(INCBLAS) \
			 -I$(ASCOT_DIR)/src/iers \
			 -I$(ASCOT_DIR)/src/ivg \
			 -I$(ASCOT_DIR)/src/ivg/lpSked \
			 -I$(ASCOT_DIR)/include \
			 -I$(ASCOT_DIR)/src/projectpluto/include \
			 -I$(ASCOT_DIR)/src/sofa/include \
			 -I$(ASCOT_DIR)/src/qt_plot \
			 -I$(ASCOT_DIR)/src/analysistools \
			 -I$(GUROBI_HOME)/include \
			 $(INCMPI)

IGG_VLBI_LIB = -L$(ASCOT_DIR)/lib -lIvgUtils -lqt_plot -lAnalysisTools $(ASCOT_DIR)/src/iers/libIers2010.a -L/usr/lib/x86_64-linux-gnu/lapack 

LIBS =						\
	-lcurl -lconfig++ \
	$(IGG_VLBI_LIB) \
	-L$(ASCOT_DIR)/src/projectpluto/lib -ljpleph -lsat \
	-L$(ASCOT_DIR)/src/sofa/lib/ -lsofa_c \
	$(LIBBLAS) \
        -lmatio -lhdf5 -L/usr/include/mysql -lmysqlclient \
	$(LIBQT) \
	-lgmp \
	-lnetcdf_c++ -lnetcdf -lfftw3 `gsl-config --libs` -lboost_regex \
	-L$(GUROBI_HOME)/lib -lgurobi_c++ -lgurobi91


#IGG_VLBI_LIB = $(LIB_DIR)/libIvgUtils.a  $(ASCOT_DIR)/lib/libIers2010.a
#
#LIBS =						\
#	-lgsl -lcurl -lconfig++ \
#	$(IGG_VLBI_LIB) \
#	-L$(ASCOT_DIR)/src/projectpluto/lib -ljpleph \
#	-L$(ASCOT_DIR)/src/sofa/lib/ -lsofa_c \
#	$(LIBBLAS) \
#	$(LIBLAPACK) -lmatio -lhdf5 -L/usr/include/mysql -lmysqlclient


all: indep_main 

indep_main: indep_main.cpp
	$(CXX) $(CXXOPTS) -o $(BIN_DIR)/ivg_indep indep_main.cpp $(LIBS) -fopenmp $(INC_ALL)
	
clean:
	rm -f $(BIN_DIR)/ivg_indep








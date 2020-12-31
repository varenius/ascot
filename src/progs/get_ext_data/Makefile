CXX=g++
CXXOPTS=-O3 -std=c++0x -fPIC

LIB_DIR = $(ASCOT_DIR)/lib
BIN_DIR = $(ASCOT_DIR)/bin

INCBLAS  = -I$(OB_DIR)/include
LIBBLAS  = -L$(OB_DIR)/lib/ -llapack -lopenblas  -lpthread -lgfortran

INCQT = -I.  -I/usr/include
LIBQT = -lQt5Gui -L/usr/lib64 -L/usr/X11R6/lib -lQt5Core -L/usr/lib/x86_64-linux-gnu/hdf5/serial -L/data1/tobson/soft/ASCOT/gurobi800/linux64/lib -L/data1/tobson/soft/ASCOT/qt5inst/qtbase/lib -lQt5Widgets -lQt5PrintSupport -lboost_program_options

#INC_ALL = -I./ $(INCQT) \
#	$(INCBLAS) \
#	-I$(ASCOT_DIR)/src/ivg \
#	-I$(ASCOT_DIR)/include
#
#LIBS =	-lgsl -lcurl -lconfig++ \
#	$(LIB_DIR)/libIvgUtils.a \
#	-L$(ASCOT_DIR)/src/projectpluto/lib -ljpleph \
#	-L$(ASCOT_DIR)/src/sofa/lib/ -lsofa_c \
#	$(LIBBLAS) \
#	$(LIBLAPACK) -lmatio -lhdf5 -L/usr/include/mysql -lmysqlclient
	
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
#			 -I$(ALGLIB_DIR)

IGG_VLBI_LIB = -L$(ASCOT_DIR)/lib -lIvgUtils -lqt_plot -lAnalysisTools $(ASCOT_DIR)/src/iers/libIers2010.a -L/usr/lib/x86_64-linux-gnu/lapack 

LIBS =						\
	-lcurl \
	$(IGG_VLBI_LIB) \
	-L$(ASCOT_DIR)/src/projectpluto/lib -ljpleph -lsat\
	-L$(ASCOT_DIR)/src/sofa/lib/ -lsofa_c \
	$(LIBBLAS) \
        -lmatio -lhdf5 -L/usr/include/mysql -lmysqlclient \
	$(LIBQT) \
	-lgmp -lCGAL  \
	-lnetcdf_c++ -lnetcdf -lfftw3 `gsl-config --libs`  -lboost_regex -lconfig++

all: get_ext_data

get_ext_data: main.cpp
	$(CXX) $(CXXOPTS) -o $(BIN_DIR)/get_ext_data main.cpp $(LIBS) $(INC_ALL) 
	
clean:
	rm -f $(BIN_DIR)/get_ext_data

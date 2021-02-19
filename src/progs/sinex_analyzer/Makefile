CXX=g++
CXXOPTS=-O3 -std=c++14 -fPIC -DBOOST_MATH_DISABLE_FLOAT128

DEFINES=-DQT_NO_DEBUG -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED

LIB_DIR = $(ASCOT_DIR)/lib
BIN_DIR = $(ASCOT_DIR)/bin

INCBLAS  = -I$(OB_DIR)/include
LIBBLAS  = -L$(OB_DIR)/lib/ -llapack -lopenblas  -lpthread -lgfortran

INCQT = -I. -I/usr/include -I/usr/include/x86_64-linux-gnu/qt5/QtWidgets/ -I/usr/include/x86_64-linux-gnu/qt5/QtCore -I/usr/include/x86_64-linux-gnu/qt5/QtGui -I/usr/include/x86_64-linux-gnu/qt5/QtPrintSupport
LIBQT = -lQt5Gui -L/usr/lib64 -L/usr/X11R6/lib -lQt5Core -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lQt5Widgets -lQt5PrintSupport -lboost_program_options 

INC_ALL = -I./ $(INCQT) \
			 $(INCBLAS) \
			 -I$(ASCOT_DIR)/src/ivg \
			 -I$(ASCOT_DIR)/include \
			 -I$(ASCOT_DIR)/src/analysistools \
			 -I$(ASCOT_DIR)/src/projectpluto/include \
			 -I$(ASCOT_DIR)/src/sofa/include \
			 -I$(ASCOT_DIR)/src/qt_plot

IGG_VLBI_LIB = -L$(ASCOT_DIR)/src/iers -lIers2010 -L$(ASCOT_DIR)/lib -lIvgUtils -lqt_plot -lAnalysisTools  $(ASCOT_DIR)/src/iers/libIers2010.a

LIBS =						\
	-lgsl -lcurl -lconfig++ \
	$(IGG_VLBI_LIB) \
	-L$(ASCOT_DIR)/src/projectpluto/lib -ljpleph -lsat\
	-L$(ASCOT_DIR)/src/sofa/lib/ -lsofa_c \
	$(LIBBLAS) \
	$(LIBLAPACK) -lmatio -lhdf5 -L/usr/include/mysql -lmysqlclient \
	$(LIBQT) \
	-lgmp -lgsl \
	-lnetcdf_c++ -lnetcdf -lfftw3 -lboost_regex

all: sinex_analyzer

sinex_analyzer: sinex_analyzer.cpp
	$(CXX) $(CXXOPTS) -o $(BIN_DIR)/sinex_analyzer sinex_analyzer.cpp $(LIBS) $(INC_ALL) 
	
clean:
	rm -f $(BIN_DIR)/sinex_analyzer








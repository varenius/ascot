CXX=g++
CXXOPTS=-O3 -std=c++14 -fPIC -DBOOST_MATH_DISABLE_FLOAT128

DEFINES=-DQT_NO_DEBUG -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED

LIB_DIR = $(ASCOT_DIR)/lib
BIN_DIR = $(ASCOT_DIR)/bin

INCBLAS  = -I$(OB_DIR)/include
LIBBLAS  = -L$(OB_DIR)/lib/ -llapack -lopenblas  -lpthread -lgfortran

INCQT = -I. -I/usr/include -I/usr/include/x86_64-linux-gnu/qt5/QtGui -I/usr/include/x86_64-linux-gnu/qt5/QtCore -I/usr/include/x86_64-linux-gnu/qt5/QtWidgets -I/usr/include/x86_64-linux-gnu/qt5/QtPrintSupport
LIBQT = -lQt5Gui -L/usr/lib64 -L/usr/X11R6/lib -lQt5Core -L/usr/lib/x86_64-linux-gnu/hdf5/serial -L/usr/include/x86_64-linux-gnu/qt5/ -lQt5Widgets -lQt5PrintSupport -lboost_program_options

INC_ALL = -I./ $(INCQT) \
			 $(INCBLAS) \
			 -I$(ASCOT_DIR)/src/ivg \
			 -I$(ASCOT_DIR)/src/analysistools \
			 -I$(ASCOT_DIR)/include \
			 -I$(ASCOT_DIR)/src/projectpluto/include \
			 -I$(ASCOT_DIR)/src/sofa/include \
			 -I$(ASCOT_DIR)/src/qt_plot

#IGG_VLBI_LIB = -L$(ASCOT_DIR)/lib -lIers2010 -lIvgUtils -lAnalysisTools -lqt_plot $(ASCOT_DIR)/src/iers/libIers2010.a
IGG_VLBI_LIB = $(LIB_DIR)/libIvgUtils.a $(LIB_DIR)/libqt_plot.a $(LIB_DIR)/libAnalysisTools.a $(ASCOT_DIR)/src/iers/libIers2010.a
#IGG_VLBI_LIB = $(LIB_DIR)/libIvgUtils.a  $(ASCOT_DIR)/lib/libIers2010.a

LIBS =						\
	-lgsl -lcurl -lconfig++ \
	$(IGG_VLBI_LIB) \
	-L$(ASCOT_DIR)/src/projectpluto/lib -ljpleph -lsat \
	-L$(ASCOT_DIR)/src/sofa/lib/ -lsofa_c \
	$(LIBBLAS) \
	$(LIBLAPACK) -lmatio -lhdf5 -L/usr/include/mysql -lmysqlclient \
	$(LIBQT) \
	-lgmp -lCGAL   -lmatio -lgsl \
	-lnetcdf_c++ -lnetcdf -lfftw3 -lboost_regex


all: analysis_main analysis_obs_main analysis_tropo_main analysis_clock_main analysis_eop_main ivg_comp_tropo

analysis_main: analysis_main.cpp
	$(CXX) $(CXXOPTS) -o $(BIN_DIR)/ivg_atools analysis_main.cpp $(LIBS) $(INC_ALL) 

analysis_obs_main: analysis_obs_main.cpp
	$(CXX) $(CXXOPTS) -o $(BIN_DIR)/ivg_atools_obs analysis_obs_main.cpp $(LIBS) $(INC_ALL) 	

analysis_tropo_main: analysis_tropo_main.cpp
	$(CXX) $(CXXOPTS) -o $(BIN_DIR)/ivg_atools_tropo analysis_tropo_main.cpp $(LIBS) $(INC_ALL) 

analysis_clock_main: analysis_clock_main.cpp
	$(CXX) $(CXXOPTS) -o $(BIN_DIR)/ivg_atools_clock analysis_clock_main.cpp $(LIBS) $(INC_ALL)
	
analysis_eop_main: analysis_eop_main.cpp
	$(CXX) $(CXXOPTS) -o $(BIN_DIR)/ivg_atools_eop analysis_eop_main.cpp $(LIBS) $(INC_ALL) 
	
ivg_comp_tropo: tropo_compare_main.cpp
	$(CXX) $(CXXOPTS) -o $(BIN_DIR)/ivg_comp_tropo tropo_compare_main.cpp $(LIBS) $(INC_ALL)
	
zwd_sinex_analyser: zwd_sinex_analyser.cpp
	$(CXX) $(CXXOPTS) -o $(BIN_DIR)/zwd_sinex_analyser zwd_sinex_analyser.cpp $(LIBS) $(INC_ALL) 
	
clean:
	rm -f $(BIN_DIR)/ivg_atools $(BIN_DIR)/ivg_atools_eop $(BIN_DIR)/ivg_atools_obs $(BIN_DIR)/ivg_atools_tropo $(BIN_DIR)/ivg_comp_tropo $(BIN_DIR)/zwd_sinex_analyser








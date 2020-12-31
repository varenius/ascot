all: analysisTools

DEFINES=
#DEFINES= -DDEBUG_ANALYSIS=1

CXX=g++
CXXOPTS=-O3 -std=c++0x $(DEFINES) -fPIC

LIBDIR = $(ASCOT_DIR)/lib
IVGUTILS = .

INCBLAS  = -I$(OB_DIR)/include
LIBBLAS  = -L$(OB_DIR)/lib/ -llapack -lopenblas  -lpthread -lgfortran


#INCQT = -I/usr/share/qt4/mkspecs/default -I. -I/usr/include/QtCore -I/usr/include/QtGui -I/usr/include
#LIBQT = -lQtGui -L/usr/lib64 -L/usr/X11R6/lib -lQtCore -L/usr/local/qt_4.x.y/lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/lapack
INCQT = -I/usr/include -I/usr/include/x86_64-linux-gnu/qt5/ -I/usr/include/x86_64-linux-gnu/qt5/QtWidgets -I/usr/include/x86_64-linux-gnu/qt5/QtCore/ -I/usr/include/x86_64-linux-gnu/qt5/QtGui -I/usr/include/x86_64-linux-gnu/qt5/QtPrintSupport
LIBQT = -L/usr/include/x86_64-linux-gnu

INC_ALL = -I./ $(INCQT) \
                         $(INCBLAS) \
                         -I$(ASCOT_DIR)/src/ivg \
			 -I$(ASCOT_DIR)/include \
			 -I$(ASCOT_DIR)/src/qt_plot \
			 -I$(ASCOT_DIR)/src/projectpluto/include \
			 -I$(ASCOT_DIR)/src/sofa/include

LIBS =                                          \
        -lgsl -lcurl -lconfig++ \
        -L$(ASCOT_DIR)/lib -lIers2010 -lIvgUtils -lqt_plot \
	-L$(ASCOT_DIR)/src/iers/libIers2010.a \
        $(LIBBLAS) \
        -lmatio -lhdf5 -L/usr/include/mysql -lmysqlclient \
	$(LIBQT)


OBJS = source_analysis.o station_analysis.o obs_analysis.o fit.o tsa.o transformation.o

source_analysis.o: source_analysis.cpp source_analysis.h
	$(CXX) $(CXXOPTS) $(OPTIMIZE) -c source_analysis.cpp $(LIBS) $(INC_ALL)

station_analysis.o: station_analysis.cpp station_analysis.h
	$(CXX) $(CXXOPTS) $(OPTIMIZE) -c station_analysis.cpp $(LIBS) $(INC_ALL)

obs_analysis.o: obs_analysis.cpp obs_analysis.h
	$(CXX) $(CXXOPTS) $(OPTIMIZE) -c obs_analysis.cpp $(LIBS) $(INC_ALL)	
	
fit.o: fit.cpp fit.h
	$(CXX) $(CXXOPTS) $(OPTIMIZE) -c fit.cpp $(LIBS) $(INC_ALL)
	
tsa.o: tsa.cpp tsa.h
	$(CXX) $(CXXOPTS) $(OPTIMIZE) -c tsa.cpp $(LIBS) $(INC_ALL)
	
transformation.o: transformation.cpp transformation.h
	$(CXX) $(CXXOPTS) $(OPTIMIZE) -c transformation.cpp $(LIBS) $(INC_ALL)

analysisTools: $(OBJS)
	rm -f $(LIBDIR)/libAnalysisTools.a
	ar rs $(LIBDIR)/libAnalysisTools.a $(OBJS)
	
clean: 
	rm -f $(OBJS)
	rm -f $(LIBDIR)/libAnalysisTools.a

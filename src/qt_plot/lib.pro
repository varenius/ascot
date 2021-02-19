# ===========================================================================
# Project File for QT-plot lib
# ===========================================================================

LIBS         = -L$(OB_DIR)/lib/ -llapack

LIBS        +=  $(ASCOT_DIR)/lib/libIvgUtils.a \
                $(ASCOT_DIR)/lib/libAnalysisTools.a \
                -lgsl -lcurl -lconfig++ \
                -L$(OB_DIR)/lib/ -llapack -lopenblas -lpthread -lgfortran \
                -lmatio -lhdf5 \
                -L/usr/include/mysql -lmysqlclient \
                -lnetcdf_c++ -lnetcdf -lfftw3 \
                -L$(GUROBI_HOME)/lib -lgurobi_c++ -lgurobi80

INCLUDEPATH +=  $(ASCOT_DIR)/src/ivg \
                $(ASCOT_DIR)/src/ivg/lpSked \
                $(ASCOT_DIR)/src/analysistools \
                $(OB_DIR)/include \
                $(ASCOT_DIR)/src/projectpluto/include \
                $(ASCOT_DIR)/src/sofa/include \
                $(ASCOT_DIR)/include  \
                $(GUROBI_HOME)/include

QT       += core gui 
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

QMAKE_CXXFLAGS += -std=c++14 -DBOOST_MATH_DISABLE_FLOAT128

TEMPLATE = lib

CONFIG += qt
CONFIG += warn_off
CONFIG += staticlib

TARGET = $(ASCOT_DIR)/lib/qt_plot

HEADERS += qcustomplot.h exportdialog.h statistics.h graphdialog.h plot.h analyzer.h projection.h \
           eopseries.h ascot.h ascotizer.h movfise.h lp_sked_gui.h lp_sked_station.h lp_sked_worker.h

SOURCES += qcustomplot.cpp exportdialog.cpp statistics.cpp graphdialog.cpp plot.cpp analyzer.cpp \
           projection.cpp eopseries.cpp  ascot.cpp ascotizer.cpp movfise.cpp lp_sked_gui.cpp lp_sked_station.cpp lp_sked_worker.cpp 

FORMS   += graphstyle.ui export.ui statistics.ui axisrange.ui plot.ui analyzer.ui projection.ui \
           eopseries.ui highlight.ui ascot.ui movfise.ui statistics_info.ui lp_sked_station.ui

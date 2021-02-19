CXX=g++
CXXOPTS=-O3 -std=c++14 -fPIC -DBOOST_MATH_DISABLE_FLOAT128

LIB_DIR = $(ASCOT_DIR)/lib
BIN_DIR = $(ASCOT_DIR)/bin

INCBLAS  = -I$(OB_DIR)/include
LIBBLAS  = -L$(OB_DIR)/lib/ -llapack -lopenblas  -lpthread -lgfortran

INC_ALL = -I./ \
	$(INCBLAS) \
	-I$(ASCOT_DIR)/src/ivg \
	-I$(ASCOT_DIR)/include \
	-I$(ASCOT_DIR)/src/projectpluto/include \
	-I$(ASCOT_DIR)/src/sofa/include

LIBS =						\
	-lgsl -lcurl -lconfig++ \
	-L$(ASCOT_DIR)/lib -lIvgUtils -lqt_plot $(ASCOT_DIR)/src/iers/libIers2010.a \
	-L$(ASCOT_DIR)/src/projectpluto/lib -ljpleph  -lsat \
	-L$(ASCOT_DIR)/src/sofa/lib/ -lsofa_c \
	$(LIBBLAS) \
	-lQt5Gui -lQt5Core -lQt5Widgets -lQt5PrintSupport \
        -L/usr/lib/x86_64-linux-gnu/hdf5/serial \
	-lgmp -lgsl  \
	-lmatio -lhdf5 -L/usr/include/mysql -lmysqlclient

all: arc_file_generator

arc_file_generator: main.cpp
	$(CXX) $(CXXOPTS) -o $(BIN_DIR)/arc_file_generator main.cpp $(LIBS) $(INC_ALL) 
	
clean:
	rm -f $(BIN_DIR)/arc_file_generator








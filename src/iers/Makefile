# ============================================================================
# makefile for fortran subroutines according IERS Conventions
#    http://maia.usno.navy.mil/conv2010/software.html
#
# build IERS-library (libiers.a)
#    make all
#
# and C++ test routines
#    make test
#
# 2014-12-03 - TA
# 2015-07-07 - TA - added IERS_CMP_2015
# 2015-10-01 - TA - added parametric
# ============================================================================
#
# you might want to change some paths
TEST_DIR = ./test
LIBDIR = ./
INCDIR = $(ASCOT_DIR)/include/

# ============================================================================
# build objects and library
IERS_OBJS = GCONV2.o           \
            PMSDNUT2.o         \
            FUNDARG.o          \
            UTLIBR.o           \
            FCNNUT.o           \
            DEHANTTIDEINEL.o   \
            ARG2.o             \
				HARDISP.o          \
            CNMTX.o            \
				ORTHO_EOP.o        \
				RG_ZONT2.o         \
				GMF.o              \
				GPT.o              \
				GPT2.o             \
                                gpt3_1.o           \
                                vmf3_ht.o          \
				VMF1.o             \
                                vmf3.o             \
				VMF1_HT.o          \
				APG.o              \
            HF2002_IERS.o      \
				IERS_CMP_2015.o    \
				parametric.o \
            calc_hfeop.o

LIBS = libIers2010.a -lgfortran

all: HARDISP XHF2002_IERS $(IERS_OBJS) iersLib

GCONV2.o: GCONV2.F
	gfortran -c GCONV2.F

PMSDNUT2.o: PMSDNUT2.F
	gfortran -c PMSDNUT2.F

FUNDARG.o: FUNDARG.F
	gfortran -c FUNDARG.F
	
UTLIBR.o: UTLIBR.F
	gfortran -c UTLIBR.F
	
FCNNUT.o: FCNNUT.F
	gfortran -c FCNNUT.F
	
ARG2.o: ARG2.F
	gfortran -c ARG2.F

DEHANTTIDEINEL.o: DEHANTTIDEINEL.F
	$(MAKE) -C dehanttideinel
	gfortran -c DEHANTTIDEINEL.F

HARDISP.o: HARDISP_ivg.F
	$(MAKE) -C hardisp
	gfortran -o HARDISP.o -c HARDISP_ivg.F
	
ORTHO_EOP.o: ORTHO_EOP.F
	gfortran -c ORTHO_EOP.F
	
CNMTX.o: CNMTX.F
	gfortran -c CNMTX.F
	
RG_ZONT2.o: RG_ZONT2.F
	gfortran -c RG_ZONT2.F
	
GMF.o: GMF.F
	gfortran -c GMF.F

GPT.o: GPT.F
	gfortran -c GPT.F

GPT2.o: GPT2.F
	gfortran -c GPT2.F

gpt3_1.o: gpt3_1.f90
	gfortran -c -ffree-line-length-512 gpt3_1.f90

VMF1.o: VMF1.F
	gfortran -c VMF1.F

VMF1_HT.o: VMF1_HT.F
	gfortran -c VMF1_HT.F

vmf3_ht.o: vmf3_ht.f90
	gfortran -c -ffree-line-length-512 vmf3_ht.f90

vmf3.o: vmf3.f90
	gfortran -c -ffree-line-length-512 vmf3.f90

APG.o: APG.F
	gfortran -c APG.F

HF2002_IERS.o: HF2002_IERS.F
	gfortran -c HF2002_IERS.F

IERS_CMP_2015.o: IERS_CMP_2015.F
	gfortran -c IERS_CMP_2015.F

parametric.o: parametric.f
	gfortran -c parametric.f

calc_hfeop.o: jd_to_date_1.f  hfeop_xyu_mod.f calc_hfeop.f
	gfortran -c -static -g -ffree-form jd_to_date_1.f
	gfortran -c -static -g -ffree-form hfeop_xyu_mod.f
	gfortran -c -static -g -ffree-form calc_hfeop.f

# ============================================================================
# build Hardisp program
HARDISP: HARDISP.F
	$(MAKE) -C hardisp
	gfortran -o HARDISP HARDISP.F libiers-hardisp.a

XHF2002_IERS: HF2002_IERS.o XHF2002_IERS.F
	gfortran -o XHF2002_IERS XHF2002_IERS.F HF2002_IERS.o

# ============================================================================
# link library
iersLib: $(OBJS)
	rm -f $(LIBDIR)/libIers2010.a
	ar -rcsT $(LIBDIR)/libIers2010.a libiers-dehant.a libiers-hardisp.a *.o 
	cp iers_wrapper.h $(INCDIR)

# ============================================================================
# build test functions/objects
TEST_CPP = $(TEST_DIR)/gconv2_test.cpp         \
           $(TEST_DIR)/pmsdnut2_test.cpp       \
           $(TEST_DIR)/utlibr_test.cpp         \
           $(TEST_DIR)/fcnnut_test.cpp         \
           $(TEST_DIR)/dehanttideinel_test.cpp \
           $(TEST_DIR)/hardisp_test.cpp            \
           $(TEST_DIR)/ortho_eop_test.cpp          \
           $(TEST_DIR)/rg_zont2_test.cpp           \
           $(TEST_DIR)/gmf_test.cpp                \
           $(TEST_DIR)/gpt_test.cpp                \
           $(TEST_DIR)/gpt2_test.cpp               \
           $(TEST_DIR)/apg_test.cpp                \
           $(TEST_DIR)/hf2002_test.cpp   

test:  $(TEST_CPP)
	g++ -std=c++11 -o $(TEST_DIR)/gconv2_test $(TEST_DIR)/gconv2_test.cpp $(LIBS) -I$(INCDIR)
	g++ -std=c++11 -o $(TEST_DIR)/pmsdnut2_test $(TEST_DIR)/pmsdnut2_test.cpp $(LIBS) -I$(INCDIR)
	g++ -std=c++11 -o $(TEST_DIR)/utlibr_test   $(TEST_DIR)/utlibr_test.cpp $(LIBS) -I$(INCDIR)
	g++ -std=c++11 -o $(TEST_DIR)/fcnnut_test   $(TEST_DIR)/fcnnut_test.cpp $(LIBS) -I$(INCDIR)
	g++ -std=c++11 -o $(TEST_DIR)/dehanttideinel_test $(TEST_DIR)/dehanttideinel_test.cpp \
                     $(LIBS) -I$(INCDIR)
	g++ -std=c++11 -o $(TEST_DIR)/hardisp_test  $(TEST_DIR)/hardisp_test.cpp -std=c++11 $(LIBS) -I$(INCDIR)
	g++ -std=c++11 -o $(TEST_DIR)/ortho_eop_test   $(TEST_DIR)/ortho_eop_test.cpp $(LIBS)  -I$(INCDIR)
	g++ -std=c++11 -o $(TEST_DIR)/rg_zont2_test   $(TEST_DIR)/rg_zont2_test.cpp $(LIBS) -I$(INCDIR)
	g++ -std=c++11 -o $(TEST_DIR)/gmf_test   $(TEST_DIR)/gmf_test.cpp $(LIBS) -I$(INCDIR)
	g++ -std=c++11 -o $(TEST_DIR)/gpt_test   $(TEST_DIR)/gpt_test.cpp $(LIBS) -I$(INCDIR)
	g++ -std=c++11 -o $(TEST_DIR)/gpt2_test   $(TEST_DIR)/gpt2_test.cpp $(LIBS) -I$(INCDIR)
	g++ -std=c++11 -o $(TEST_DIR)/apg_test   $(TEST_DIR)/apg_test.cpp $(LIBS) -I$(INCDIR)
	g++ -std=c++11 -o $(TEST_DIR)/hf2002_test   $(TEST_DIR)/hf2002_test.cpp -std=c++11 $(LIBS) -I$(INCDIR)

# ============================================================================
clean:
	rm -f $(IERS_OBJS) HARDISP XHF2002_IERS $(LIBDIR)/libIers2010.a
	$(MAKE) clean -C dehanttideinel
	$(MAKE) clean -C hardisp

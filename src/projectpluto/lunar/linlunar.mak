# GNU MAKE Makefile for 'lunar' basic astronomical functions library
#
# Usage: make -f [path\]linlunar.mak [CLANG=Y] [XCOMPILE=Y] [MSWIN=Y] [tgt]
#
# where tgt can be any of:
# [all|astcheck|astephem|calendar... clean]
#
#	'XCOMPILE' = cross-compile for Windows,  using MinGW,  on a Linux box
#	'MSWIN' = compile for Windows,  using MinGW,  on a Windows machine
#	'CLANG' = use clang instead of GCC;  Linux only
# None of these: compile using g++ on Linux,  for Linux

CC=g++
LIBSADDED=
EXE=

ifdef CLANG
	CC=clang
endif

RM=rm -f

ifdef MSWIN
	EXE=.exe
else
	LIBSADDED=-lm -lrt
endif

ifdef XCOMPILE
	CC=x86_64-w64-mingw32-g++
	EXE=.exe
	LIBSADDED=
endif

all: astcheck$(EXE) astephem$(EXE) calendar$(EXE) colors$(EXE) colors2$(EXE) \
	cosptest$(EXE) dist$(EXE) easter$(EXE) get_test$(EXE) htc20b$(EXE) jd$(EXE) \
 jevent$(EXE) jpl2b32$(EXE) jsattest$(EXE) lun_test$(EXE) marstime$(EXE) oblitest$(EXE) persian$(EXE)  \
 phases$(EXE) ps_1996$(EXE) ssattest$(EXE) tables$(EXE) \
	test_ref$(EXE) testprec$(EXE) uranus1$(EXE) utc_test$(EXE)

CFLAGS=-Wall -O3

.cpp.o:
	$(CC) $(CFLAGS) -c $<

OBJS= alt_az.o astfuncs.o big_vsop.o classel.o cospar.o date.o delta_t.o \
	de_plan.o dist_pa.o eart2000.o elp82dat.o getplane.o get_time.o \
	jsats.o lunar2.o miscell.o nutation.o obliquit.o pluto.o precess.o \
	showelem.o ssats.o triton.o vsopson.o

lunar.a: $(OBJS)
	$(RM) lunar.a
	ar rv lunar.a $(OBJS)

clean:
	$(RM) $(OBJS)
	$(RM) astcheck.o astephem.o calendar.o cosptest.o get_test.o gust86.o
	$(RM) htc20b.o jd.o jevent.o jpl2b32.o jsattest.o lun_test.o
	$(RM) lun_tran.o mpcorb.o oblitest.o obliqui2.o persian.o phases.o
	$(RM) ps_1996.o refract.o refract4.o riseset3.o solseqn.o spline.o
	$(RM) ssattest.o tables.o test_ref.o testprec.o uranus1.o utc_test.o
	$(RM) astcheck$(EXE) astephem$(EXE) calendar$(EXE) colors$(EXE)
	$(RM) colors2$(EXE) cosptest$(EXE) dist$(EXE) easter$(EXE) get_test$(EXE)
	$(RM) htc20b$(EXE) jd$(EXE) jevent$(EXE) jpl2b32$(EXE) jsattest$(EXE)
	$(RM) lun_test$(EXE) marstime$(EXE) oblitest$(EXE) persian$(EXE) phases$(EXE)
	$(RM) ps_1996$(EXE) relativi$(EXE) ssattest$(EXE) tables$(EXE)
	$(RM) test_ref$(EXE) testprec$(EXE) uranus1$(EXE) utc_test$(EXE) lunar.a

astcheck$(EXE): astcheck.o mpcorb.o lunar.a
	$(CC) $(CFLAGS) -o astcheck$(EXE) astcheck.o mpcorb.o lunar.a $(LIBSADDED)

astephem$(EXE): astephem.o mpcorb.o lunar.a
	$(CC) $(CFLAGS) -o astephem$(EXE) astephem.o mpcorb.o lunar.a $(LIBSADDED)

calendar$(EXE): calendar.o lunar.a
	$(CC) $(CFLAGS) -o calendar$(EXE) calendar.o   lunar.a $(LIBSADDED)

colors$(EXE): colors.cpp
	$(CC) $(CFLAGS) -o colors$(EXE) colors.cpp -DSIMPLE_TEST_PROGRAM

colors2$(EXE): colors2.cpp
	$(CC) $(CFLAGS) -o colors2$(EXE) colors2.cpp -DTEST_FUNC

cosptest$(EXE): cosptest.o lunar.a
	$(CC) $(CFLAGS) -o cosptest$(EXE) cosptest.o   lunar.a $(LIBSADDED)

dist$(EXE): dist.cpp
	$(CC) $(CFLAGS) -o dist$(EXE) dist.cpp $(LIBSADDED)

easter$(EXE): easter.cpp lunar.a
	$(CC) $(CFLAGS) -o easter$(EXE) -DTEST_CODE easter.cpp lunar.a $(LIBSADDED)

get_test$(EXE): get_test.o lunar.a
	$(CC) $(CFLAGS) -o get_test$(EXE) get_test.o lunar.a $(LIBSADDED)

htc20b$(EXE): htc20b.cpp lunar.a
	$(CC) $(CFLAGS) -o htc20b$(EXE) -DTEST_MAIN htc20b.cpp lunar.a $(LIBSADDED)

jd$(EXE): jd.o lunar.a
	$(CC) $(CFLAGS) -o jd$(EXE) jd.o lunar.a $(LIBSADDED)

jevent$(EXE):	jevent.o lunar.a
	$(CC) $(CFLAGS) -o jevent$(EXE) jevent.o lunar.a $(LIBSADDED)

jpl2b32$(EXE):	jpl2b32.o
	$(CC) $(CFLAGS) -o jpl2b32$(EXE) jpl2b32.o $(LIBSADDED)

jsattest$(EXE): jsattest.o lunar.a
	$(CC) $(CFLAGS) -o jsattest$(EXE) jsattest.o lunar.a $(LIBSADDED)

lun_test$(EXE): lun_test.o lun_tran.o riseset3.o lunar.a
	$(CC)	$(CFLAGS) -o lun_test$(EXE) lun_test.o lun_tran.o riseset3.o lunar.a $(LIBSADDED)

marstime$(EXE): marstime.cpp
	$(CC) $(CFLAGS) -o marstime$(EXE) marstime.cpp -DTEST_PROGRAM $(LIBSADDED)

oblitest$(EXE): oblitest.o obliqui2.o spline.o lunar.a
	$(CC) $(CFLAGS) -o oblitest$(EXE) oblitest.o obliqui2.o spline.o lunar.a $(LIBSADDED)

persian$(EXE): persian.o solseqn.o lunar.a
	$(CC) $(CFLAGS) -o persian$(EXE) persian.o solseqn.o lunar.a $(LIBSADDED)

phases$(EXE): phases.o lunar.a
	$(CC) $(CFLAGS) -o phases$(EXE)   phases.o   lunar.a $(LIBSADDED)

ps_1996$(EXE): ps_1996.o lunar.a
	$(CC) $(CFLAGS) -o ps_1996$(EXE)   ps_1996.o   lunar.a $(LIBSADDED)

relativi$(EXE): relativi.cpp lunar.a
	$(CC) $(CFLAGS) -o relativi$(EXE) -DTEST_CODE relativi.cpp lunar.a $(LIBSADDED)

ssattest$(EXE): ssattest.o lunar.a
	$(CC) $(CFLAGS) -o ssattest$(EXE) ssattest.o lunar.a $(LIBSADDED)

tables$(EXE):                    tables.o riseset3.o lunar.a
	$(CC) $(CFLAGS) -o tables$(EXE) tables.o riseset3.o lunar.a $(LIBSADDED)

test_ref$(EXE):                    test_ref.o refract.o refract4.o
	$(CC) $(CFLAGS) -o test_ref$(EXE) test_ref.o refract.o refract4.o $(LIBSADDED)

testprec$(EXE):                    testprec.o lunar.a
	$(CC) $(CFLAGS) -o testprec$(EXE) testprec.o lunar.a $(LIBSADDED)

uranus1$(EXE): uranus1.o gust86.o
	$(CC) $(CFLAGS) -o uranus1$(EXE) uranus1.o gust86.o $(LIBSADDED)

utc_test$(EXE):                utc_test.o lunar.a
	$(CC) $(CFLAGS) -o utc_test$(EXE) utc_test.o lunar.a $(LIBSADDED)


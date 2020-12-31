# Makefile for gcc (and MinGW,  and clang)
# Note dependence of 'sub_eph' on the 'lunar' library.  This is available
# at http://www.projectpluto.com/source.htm .
# Note also that 'sub_eph' is the last file to be built.  You can make
# everything else and ignore the missing dependency error if you aren't
# interested in building 'sub_eph'.

ifdef CLANG
CC=clang
LIB=-lm -L../lib -llunar
else
CC=g++
LIB= -L../lib -llunar
endif

ifdef MSWIN
EXE=.exe
else
EXE=
endif

all: asc2eph$(EXE) dump_eph$(EXE) eph2asc$(EXE) merge_de$(EXE) testeph$(EXE) sub_eph$(EXE)

OBJS=f_strtod.o  jpleph.o

CFLAGS=-c -Wall -O3
EPHLIB=../lib/libjpleph.a

.cpp.o:
	$(CC) $(CFLAGS) $<

$(EPHLIB): $(OBJS)
	rm -f $(EPHLIB)
	ar rv $(EPHLIB) $(OBJS)
	cp jpleph.h ../include

asc2eph$(EXE):    asc2eph.o f_strtod.o $(EPHLIB)
	$(CC) -o asc2eph asc2eph.o f_strtod.o $(LIB)  $(EPHLIB)

eph2asc$(EXE):    eph2asc.o jpleph.o $(EPHLIB)
	$(CC) -o eph2asc eph2asc.o jpleph.o $(LIB) $(EPHLIB)

dump_eph$(EXE):    dump_eph.o jpleph.o $(EPHLIB)
	$(CC) -o dump_eph dump_eph.o jpleph.o $(LIB) $(EPHLIB)

merge_de$(EXE):    merge_de.o jpleph.o $(EPHLIB)
	$(CC) -o merge_de merge_de.o jpleph.o $(LIB) $(EPHLIB)

sub_eph$(EXE):    sub_eph.o jpleph.o  $(EPHLIB)
	$(CC) -o sub_eph sub_eph.o jpleph.o $(LIB)  $(EPHLIB) -I../include -I./

sub_eph.o: sub_eph.cpp $(EPHLIB)
	$(CC) $(CFLAGS) -DTEST_MAIN sub_eph.cpp -I../include $(EPHLIB)

testeph$(EXE):    testeph.o jpleph.o $(EPHLIB)
	$(CC) -o testeph testeph.o jpleph.o $(LIB) -I../include -I./ $(EPHLIB)

clean:
	$(RM) *.o
	$(RM) asc2eph$(EXE)
	$(RM) dump_eph$(EXE)
	$(RM) eph2asc$(EXE)
	$(RM) merge_de$(EXE)
	$(RM) sub_eph$(EXE)
	$(RM) testeph$(EXE)  
	$(RM) $(EPHLIB)
	$(RM) ../include/jpleph.h

all: astcheck astephem calendar easter get_test \
 htc20b jd jsattest lun_test oblitest persian  \
 phases testprec tables \
	ssattest uranus1

CC=g++

CFLAGS=-Wall -O3

LUNLIB=../lib/liblunar.a

.cpp.o:
	$(CC) $(CFLAGS) -c $<

OBJS= alt_az.o astfuncs.o big_vsop.o classel.o cospar.o date.o delta_t.o \
	de_plan.o dist_pa.o eart2000.o elp82dat.o getplane.o get_time.o \
	jsats.o lunar2.o miscell.o nutation.o obliquit.o pluto.o precess.o \
	showelem.o ssats.o triton.o vsopson.o

$(LUNLIB): $(OBJS)
	rm -f $(LUNLIB)
	ar rv $(LUNLIB) $(OBJS)
	cp watdefs.h date.h ../include/

astcheck:  astcheck.o mpcorb.o $(LUNLIB)
	$(CC) $(CFLAGS) -o astcheck astcheck.o mpcorb.o $(LUNLIB)

astephem:  astephem.o mpcorb.o $(LUNLIB)
	$(CC) $(CFLAGS) -o astephem astephem.o mpcorb.o $(LUNLIB)

calendar: calendar.o $(LUNLIB)
	$(CC) $(CFLAGS) -o calendar   calendar.o   $(LUNLIB)

easter: easter.cpp $(LUNLIB)
	$(CC) $(CFLAGS) -o easter -DTEST_CODE easter.cpp $(LUNLIB)

get_test: get_test.o $(LUNLIB)
	$(CC) $(CFLAGS) -o get_test get_test.o $(LUNLIB)

htc20b: htc20b.cpp $(LUNLIB)
	$(CC) $(CFLAGS) -o htc20b -DTEST_MAIN htc20b.cpp $(LUNLIB)

jd: jd.o $(LUNLIB)
	$(CC) $(CFLAGS) -o jd jd.o $(LUNLIB)

jsattest: jsattest.o $(LUNLIB)
	$(CC) $(CFLAGS) -o jsattest jsattest.o $(LUNLIB)

lun_test:                lun_test.o lun_tran.o riseset3.o $(LUNLIB)
	$(CC)	$(CFLAGS) -o lun_test lun_test.o lun_tran.o riseset3.o $(LUNLIB)

oblitest: oblitest.o obliqui2.o spline.o $(LUNLIB)
	$(CC) $(CFLAGS) -o oblitest oblitest.o obliqui2.o spline.o $(LUNLIB)

persian: persian.o solseqn.o $(LUNLIB)
	$(CC) $(CFLAGS) -o persian persian.o solseqn.o $(LUNLIB)

phases: phases.o $(LUNLIB)
	$(CC) $(CFLAGS) -o phases   phases.o   $(LUNLIB)

testprec: testprec.o obliqui2.o spline.o $(LUNLIB)
	$(CC) $(CFLAGS) -o testprec testprec.o obliqui2.o spline.o $(LUNLIB)

relativi: relativi.cpp $(LUNLIB)
	$(CC) $(CFLAGS) -o relativi -DTEST_CODE relativi.cpp $(LUNLIB)

ssattest: ssattest.o $(LUNLIB)
	$(CC) $(CFLAGS) -o ssattest ssattest.o $(LUNLIB)

tables: tables.o riseset3.o $(LUNLIB)
	$(CC) $(CFLAGS) -o tables tables.o riseset3.o $(LUNLIB)

uranus1: uranus1.o gust86.o
	$(CC) $(CFLAGS) -o uranus1 uranus1.o gust86.o

clean:
	rm -f *.o
	rm -f $(LUNLIB)
	rm -f ../include/watdefs.h
	rm -f ../include/date.h 

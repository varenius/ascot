CC=g++
EXE=
RM=rm -f

ifdef CLANG
	CC=clang
endif

ifdef MSWIN
	EXE=.exe
endif

ifdef XCOMPILE
	CC=x86_64-w64-mingw32-g++
	EXE=.exe
endif

all: get_high$(EXE) mergetle$(EXE) obs_tes2$(EXE) obs_test$(EXE) out_comp$(EXE) \
	test_sat$(EXE) test2$(EXE) sat_id$(EXE) sat_id2$(EXE) test_out$(EXE)

CFLAGS=-Wall -O3 -pedantic -Wextra

clean:
	$(RM) *.o
	$(RM) get_high$(EXE)
	$(RM) mergetle$(EXE)
	$(RM) obs_tes2$(EXE)
	$(RM) obs_test$(EXE)
	$(RM) out_comp$(EXE)
	$(RM) sat_code.a
	$(RM) sat_id$(EXE)
	$(RM) sat_id2$(EXE)
	$(RM) test2$(EXE)
	$(RM) test_out$(EXE)
	$(RM) test_sat$(EXE)
	$(RM) ../lib/libsat.a

OBJS= sgp.o sgp4.o sgp8.o sdp4.o sdp8.o deep.o basics.o get_el.o common.o

get_high$(EXE):	 get_high.o get_el.o
	$(CC) $(CFLAGS) -o get_high$(EXE) get_high.o get_el.o

mergetle$(EXE):	 mergetle.o
	$(CC) $(CFLAGS) -o mergetle$(EXE) mergetle.o

obs_tes2$(EXE):	 obs_tes2.o observe.o sat_code.a
	$(CC) $(CFLAGS) -o obs_tes2$(EXE) obs_tes2.o observe.o sat_code.a -lm

obs_test$(EXE):	 obs_test.o observe.o sat_code.a
	$(CC) $(CFLAGS) -o obs_test$(EXE) obs_test.o observe.o sat_code.a -lm

out_comp$(EXE):	 out_comp.o
	$(CC) $(CFLAGS) -o out_comp$(EXE) out_comp.o -lm

sat_code.a: $(OBJS)
	rm -f sat_code.a
	ar rv ../lib/libsat.a $(OBJS)
	cp norad.h ../include

sat_id$(EXE):	 	sat_id.o	observe.o sat_code.a
	$(CC) $(CFLAGS) -o sat_id$(EXE) sat_id.o observe.o sat_code.a -lm

sat_id2$(EXE):	 	sat_id2.o sat_id.cpp observe.o sat_code.a ../find_orb/cgi_func.cpp
	$(CC) $(CFLAGS) -o sat_id2$(EXE) -DON_LINE_VERSION sat_id2.o sat_id.cpp observe.o ../find_orb/cgi_func.cpp sat_code.a -lm

test2$(EXE):	 	test2.o sgp.o sat_code.a
	$(CC) $(CFLAGS) -o test2$(EXE) test2.o sgp.o sat_code.a -lm

test_out$(EXE):	 test_out.o tle_out.o get_el.o sgp4.o common.o
	$(CC) $(CFLAGS) -o test_out$(EXE) test_out.o tle_out.o get_el.o sgp4.o common.o -lm

test_sat$(EXE):	 test_sat.o sat_code.a
	$(CC) $(CFLAGS) -o test_sat$(EXE) test_sat.o sat_code.a -lm

.cpp.o:
	$(CC) $(CFLAGS) -c $<

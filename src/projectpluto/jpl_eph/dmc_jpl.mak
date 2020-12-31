all: testeph.exe merge_de.exe dump_eph.exe asc2eph.exe eph2asc.exe

testeph.exe:      testeph.obj jpleph.obj
   link           testeph.obj jpleph.obj

eph2asc.exe:    eph2asc.obj jpleph.obj
   link         eph2asc.obj jpleph.obj

merge_de.exe: merge_de.obj jpleph.obj
   link       merge_de.obj jpleph.obj

dump_eph.exe: dump_eph.obj jpleph.obj
   link       dump_eph.obj jpleph.obj

asc2eph.exe:  asc2eph.obj f_strtod.obj
   link       asc2eph.obj f_strtod.obj

CFLAGS=-c

.cpp.obj:
   dmc $< $(CFLAGS)

jpleph.obj:

testeph.obj:

merge_de.obj:

dump_eph.obj:

asc2eph.obj:

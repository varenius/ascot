# Note dependence of 'sub_eph' on the 'lunar' library.  This is available
# at http://www.projectpluto.com/source.htm .

all: asc2eph.exe dump_eph.exe eph2asc.exe merge_de.exe sub_eph.exe testeph.exe

eph2asc.exe:     eph2asc.obj jpleph.obj
   wcl386 -zq -k20000 eph2asc.obj jpleph.obj

dump_eph.exe:     dump_eph.obj jpleph.obj
   wcl386 -zq -k20000 dump_eph.obj jpleph.obj

testeph.exe:      testeph.obj jpleph.obj
   wcl386 -zq -k20000 testeph.obj jpleph.obj

sub_eph.exe:      sub_eph.obj jpleph.obj
   wcl386 -zq -k20000 sub_eph.obj jpleph.obj wafuncs.lib

asc2eph.exe:          asc2eph.obj f_strtod.obj
   wcl386 -zq -k20000 asc2eph.obj f_strtod.obj

merge_de.exe: merge_de.obj jpleph.obj
   wcl386 -zq -k20000 merge_de.obj jpleph.obj

CFLAGS=-W4 -Ox -j -4r -s -zq

.cpp.obj:
   wcc386 $(CFLAGS) $<

jpleph.obj:

testeph.obj:

sub_eph.obj: sub_eph.cpp
   wcc386 $(CFLAGS) -DTEST_MAIN sub_eph.cpp

dump_eph.obj:

merge_de.obj:

asc2eph.obj:


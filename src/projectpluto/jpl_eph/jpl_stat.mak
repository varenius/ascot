# Visual C/C++ makefile with static (non-DLL) linking
# Note dependence of 'sub_eph' on the 'lunar' library.  This is available
# at http://www.projectpluto.com/source.htm .

all: testeph.exe merge_de.exe dump_eph.exe asc2eph.exe sub_eph.exe

testeph.exe:    testeph.obj jpleph.obj
   link /nologo testeph.obj jpleph.obj

merge_de.exe: merge_de.cpp jpleph.obj
   cl -W3 -Ox -nologo merge_de.cpp jpleph.obj

dump_eph.exe:   dump_eph.obj jpleph.obj
   link /nologo dump_eph.obj jpleph.obj

asc2eph.exe: asc2eph.obj
   link /nologo asc2eph.obj

sub_eph.exe:    sub_eph.obj jpleph.obj lunar.lib
   link /nologo sub_eph.obj jpleph.obj lunar.lib

jpleph.obj: jpleph.cpp
   cl -W3 -Ox -GX -c -nologo jpleph.cpp

testeph.obj: testeph.cpp
   cl -W3 -Ox -GX -c -nologo testeph.cpp

sub_eph.obj: sub_eph.cpp
   cl -W3 -Ox -GX -c -nologo -DTEST_MAIN sub_eph.cpp

dump_eph.obj: dump_eph.cpp
   cl -W3 -Ox -GX -c -nologo dump_eph.cpp

asc2eph.obj: asc2eph.cpp
   cl -W3 -Ox -nologo asc2eph.cpp

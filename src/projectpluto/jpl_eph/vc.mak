# Note dependence of 'sub_eph' on the 'lunar' library.  This is available
# at http://www.projectpluto.com/source.htm .

all: asc2eph.exe dump_eph.exe eph2asc.exe merge_de.exe testeph.exe sub_eph.exe

testeph.exe:    testeph.obj jpleph.dll
   link /nologo testeph.obj jpleph.lib

merge_de.exe:   merge_de.obj jpleph.lib
   link /nologo merge_de.obj jpleph.lib

eph2asc.exe:    eph2asc.obj jpleph.lib
   link /nologo eph2asc.obj jpleph.lib

dump_eph.exe:   dump_eph.obj jpleph.lib
   link /nologo dump_eph.obj jpleph.lib

asc2eph.exe:    asc2eph.obj f_strtod.obj
   link /nologo asc2eph.obj f_strtod.obj

sub_eph.exe:    sub_eph.obj jpleph.dll lunar.lib
   link /nologo sub_eph.obj jpleph.lib lunar.lib

jpleph.lib: jpleph.obj
   del jpleph.lib
   del jpleph.dll
   link /nologo /DLL /IMPLIB:jpleph.lib /DEF:jpleph.def jpleph.obj

jpleph.obj: jpleph.cpp
   cl -W3 -Ox -GX -c -LD -nologo jpleph.cpp

testeph.obj: testeph.cpp
   cl -W3 -Ox -GX -c -nologo testeph.cpp

sub_eph.obj: sub_eph.cpp
   cl -W3 -Ox -GX -c -nologo -DTEST_MAIN sub_eph.cpp

CFLAGS=-W3 -Ox -GX -c -nologo

.cpp.obj:
   cl $(CFLAGS) $< >> err
   type err


# Makefile for MSVC
all:  test2.exe test_sat.exe obs_test.exe obs_tes2.exe sat_id.exe out_comp.exe


clean:
   del *.obj
   del *.exe

LINK=link /nologo /stack:10000

test2.exe: test2.obj sgp.obj sgp4.obj sgp8.obj sdp4.obj sdp8.obj \
          deep.obj basics.obj get_el.obj common.obj
   $(LINK) /stack:10000 test2.obj sgp.obj sgp4.obj sgp8.obj sdp4.obj sdp8.obj \
          deep.obj basics.obj get_el.obj common.obj

test_sat.exe: test_sat.obj sgp.obj sgp4.obj sgp8.obj sdp4.obj sdp8.obj \
          deep.obj basics.obj get_el.obj common.obj
   $(LINK)    test_sat.obj sgp.obj sgp4.obj sgp8.obj sdp4.obj sdp8.obj \
          deep.obj basics.obj get_el.obj common.obj

obs_test.exe: obs_test.obj sgp.obj sgp4.obj sgp8.obj sdp4.obj sdp8.obj \
          deep.obj basics.obj get_el.obj observe.obj common.obj
   $(LINK)       obs_test.obj sgp.obj sgp4.obj sgp8.obj sdp4.obj sdp8.obj \
          deep.obj basics.obj get_el.obj observe.obj common.obj

obs_tes2.exe: obs_tes2.obj sgp.obj sgp4.obj sdp4.obj \
          deep.obj basics.obj get_el.obj observe.obj common.obj
   $(LINK)       obs_tes2.obj sgp.obj sgp4.obj sdp4.obj \
          deep.obj basics.obj get_el.obj observe.obj common.obj

sat_id.exe: sat_id.obj sgp.obj sgp4.obj sdp4.obj \
          deep.obj basics.obj get_el.obj observe.obj common.obj
   $(LINK)  sat_id.obj sgp.obj sgp4.obj sdp4.obj \
          deep.obj basics.obj get_el.obj observe.obj common.obj

out_comp.exe: out_comp.obj
   $(LINK)    out_comp.obj

CFLAGS=-W3 -c -Ox -D_CRT_SECURE_NO_WARNINGS -nologo
#CFLAGS=-W3 -c -Ox -DRETAIN_PERTURBATION_VALUES_AT_EPOCH

.cpp.obj:
   cl $(CFLAGS) $<


# 2D fluid simulation makefile
#
# Author : Chris H. Rycroft (Harvard SEAS / LBL)
# Email  : chr@alum.mit.edu
# Date   : June 15th 2015

# Load the common configuration file
include ../../config.mk

iflags=-I../../tgmg -I../../shared -I../../levelset
lflags=-L../../shared -L../../levelset -L.

objs=common.o fluid_2d.o fluid_2d_io.o mgs_common.o mgs_mac.o mgs_fem.o \
	 sim_type.o object.o obj_field.o bi_interp.o
src=$(patsubst %.o,%.cc,$(objs))
execs=mg_test ftest extrema tr_unpack tr_analyze conv_test sediment

all:
	$(MAKE) -C ../../shared
	$(MAKE) -C ../../tgmg lib
	$(MAKE) -C ../../levelset
	$(MAKE) executables

executables: $(execs)

depend: $(src)
	$(cxx) $(cflags) $(iflags) -MM $(src) >Makefile.dep

-include Makefile.dep

libf2d.a: $(objs)
	rm -f libf2d.a
	ar rs libf2d.a $^

mg_test: mg_test.cc libf2d.a
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $< -lf2d -llevel++

ftest: ftest.cc libf2d.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lf2d -llevel++

sediment: sediment.cc libf2d.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lf2d -llevel++

conv_test: conv_test.cc libf2d.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lf2d -llevel++

extrema: extrema.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) -lgpmtx $(png_lflags)

tr_unpack: tr_unpack.cc common.o
	$(cxx) $(cflags) $(iflags) -o $@ $^

tr_analyze: tr_analyze.cc common.o
	$(cxx) $(cflags) $(iflags) -o $@ $^

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

clean:
	rm -rf $(execs) $(objs) libf2d.a temp_* gp_headers/temp*.gnuplot

.PHONY: clean all executables depend

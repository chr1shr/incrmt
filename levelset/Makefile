#Load the common configuration file
include ../../config.mk

#List of the common source files
objs=heap.o ex_standard.o sorting.o helper.o levelset.o
src=$(patsubst %.o,%.cc,$(objs))
execs=newton_test ex_test growth

all: $(execs)

depend:
	$(cxx) -MM $(src) >Makefile.dep

-include Makefile.dep

built/bicubic.cc: built/bicubic.pl
	cd built && perl bicubic.pl

built/bicubic2.cc: built/bicubic2.pl
	cd built && perl bicubic2.pl

built/leonardo.cc built/leonardo.hh: built/leonardo.pl
	cd built && perl leonardo.pl

newton_test: newton_test.cc liblevel++.a
	$(cxx) $(cflags) -o newton_test newton_test.cc -L. -llevel++

growth: growth.cc liblevel++.a
	$(cxx) $(cflags) -o growth growth.cc -L. -llevel++

ex_test: ex_test.cc liblevel++.a
	$(cxx) $(cflags) -o ex_test ex_test.cc -L. -llevel++

exs_test: exs_test.cc liblevel++.a ../shared/vec.hh
	$(cxx) $(cflags) -o exs_test exs_test.cc -L. -llevel++

%.o: %.cc
	$(cxx) $(cflags) -c $<

liblevel++.a: $(objs)
	rm -f liblevel++.a
	ar rs liblevel++.a $^

clean:
	rm -f $(objs) $(execs) liblevel++.a
	rm -f built/bicubic.cc built/bicubic2.cc built/leonardo.cc built/leonardo.hh

#Load the common configuration file
include ../../../config.mk

execs=quicksort quicksort2 smoothsort smoothsort2 smoothsort3

all: $(execs)

quicksort: quicksort.cc
	$(cxx) $(cflags) -o $@ $<

quicksort2: quicksort2.cc
	$(cxx) $(cflags) -o $@ $<

smoothsort: smoothsort.cc
	$(cxx) $(cflags) -o $@ $<

smoothsort2: smoothsort2.cc
	$(cxx) $(cflags) -o $@ $<

smoothsort3: smoothsort3.cc leonardo.cc
	$(cxx) $(cflags) -o $@ $<

leonardo.cc: leonardo.pl
	perl leonardo.pl

clean:
	rm -f $(execs) leonardo.cc

# This Makefile is a part of the JIMWLK solver.
# Copyright (C) 2011 Bjoern Schenke.
# JIMWLK-Solver is licenced under the GNU GPL version 2, see COPYING for details.

SHELL = /bin/sh
include ../Makefile.inc

TMPDIR = /tmp
INCLUDES = -I. 
OBJLIBS = ../lib/libsrc.a 
SOURCES = FFT.cpp Matrix.cpp Setup.cpp Init.cpp Random.cpp Group.cpp Lattice.cpp Cell.cpp Measure.cpp
OBJECTS = $(SOURCES:.cpp=.o)
CXXFLAGS = -O2 `../../bin/gsl-config --cflags` 

all: $(OBJLIBS)

../lib/libsrc.a : $(OBJECTS)
	$(AR) $(ARFLAGS) ../lib/libsrc.a $(OBJECTS)
	$(RANLIB) ../lib/libsrc.a

# general rule to make a .o out of a .cpp file
.cpp.o:
	$(CXX) $(CXXFLAGS) -Wno-write-strings -c $(INCLUDES) -o $@ $<


.PHONY: clean
clean :
	-$(RM) -f $(OBJECTS)

force_look :
	true

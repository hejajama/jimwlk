SHELL = /bin/sh

include Makefile.inc

INCLUDES   = -I./src
OBJECTS	= main.o
DIRS = src
OBJLIBS	= ./lib/libsrc.a
LIBS	=  `../bin/gsl-config --libs` -L./lib  -l:libfftw3.so.3
EXECUTABLE = main
CXXFLAGS = -O2 `../bin/gsl-config --cflags`

all: main

$(EXECUTABLE): $(OBJECTS) $(OBJLIBS) 
	$(CXX) $(CXXFLAGS) $@.o $(LIBS) -o $@ -lsrc

main.o: $(OBJLIBS) main.cpp
	$(CXX) $(CXXFLAGS) -c $(INCLUDES) main.cpp

./lib/libsrc.a : force_look
	cd src; $(MAKE) $(MFLAGS)

.PHONY: clean # explicitly declare 'clean' to be PHONY
clean :
	-$(RM) -f $(EXECUTABLE) $(OBJECTS) $(OBJLIBS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look :
	true

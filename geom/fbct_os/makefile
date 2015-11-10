# This Makefile is used under Linux

MATLABDIR ?= /Applications/MATLAB_R2011a_Student.app
CXX ?= gcc
#CXX = g++-3.3
CC ?= g++
CFLAGS = -Wall -Wconversion -O3 -O2 -fPIC -I$(MATLABDIR)/extern/include -I..

MEX = $(MATLABDIR)/bin/mex
MEX_OPTION = CC\#$(CXX) CXX\#$(CXX) CFLAGS\#"$(CFLAGS)" CXXFLAGS\#"$(CFLAGS)"
# comment the following line if you use MATLAB on a 32-bit computer
MEX_OPTION += -largeArrayDims -I$(MATLABDIR)/extern/include
MEX_EXT = $(shell $(MATLABDIR)/bin/mexext)


all:	matlab

matlab:	binary

binary: fbct_geom_mex.$(MEX_EXT)


fbct_geom_mex.$(MEX_EXT): fbct_geom_mex.c fbct_proj.o fbct_back.o 
		$(MEX) $(MEX_OPTION) fbct_geom_mex.c fbct_proj.o fbct_back.o 
								
fbct_proj.o : fbct_proj.c fbct_proj.h ct_def.h
		$(CXX) $(CFLAGS) -c fbct_proj.c fbct_proj.h ct_def.h
		
fbct_back.o : fbct_back.c fbct_back.h ct_def.h
		$(CXX) $(CFLAGS) -c fbct_back.c fbct_back.h ct_def.h
		

clean:
	rm -f *~ *.o *.obj *.gch
	@echo "clean completed"

clear:
	rm -f *~ *.o *.mex* *.obj *.gch
	@echo "clean completed"
	

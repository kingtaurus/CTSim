# This Makefile is used under Linux

MATLABDIR ?= /Applications/MATLAB_R2011a_Student.app
CXX ?= g++
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

binary: cbct_geom_mex.$(MEX_EXT)


cbct_geom_mex.$(MEX_EXT): cbct_geom_mex.c cbct_proj.o cbct_back.o 
		$(MEX) $(MEX_OPTION) cbct_geom_mex.c cbct_proj.o cbct_back.o 
								
cbct_proj.o : cbct_proj.c cbct_proj.h ct_def.h
		$(CXX) $(CFLAGS) -c cbct_proj.c cbct_proj.h ct_def.h
		
cbct_back.o : cbct_back.c cbct_back.h ct_def.h
		$(CXX) $(CFLAGS) -c cbct_back.c cbct_back.h ct_def.h
		

clean:
	rm -f *~ *.o *.obj *.gch
	@echo "clean completed"

clear:
	rm -f *~ *.o *.mex* *.obj *.gch
	@echo "clean completed"
	

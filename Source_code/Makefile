#SJL: 10/15
#The make file for the HERCULES code
#made 11/12/15 using intel/15.0.0-fasrc01


CFLAGS = -O3 -std=c++11 #-o3 -ansi -g -Wall
LIBS = #-l stdc++ -lm
#CXX = g++
CXX=icc
OBJECTS= hercules.o initialisation.o EOS_functions.o vector_operators.o planet_functions.o concentric_layer_functions.o general_functions.o itterate.o parameters_functions.o io_routines.o

# These lines define what each file depends on, and how to make it
all: HERCULESv1.0 #clean #add clean if needed

HERCULESv1.0: ${OBJECTS}
	$(CXX) $(CFLAGS) $(LIBS) ${OBJECTS} -o HERCULESv1.0

%.o: %.cc
	$(CXX) $(CFLAGS) $< -c -o $@
%: %.cc Makefile
	$(CXX) $(CFLAGS) $(LIBS) $< -o $@


clean:
	rm *.o
	echo Clean done

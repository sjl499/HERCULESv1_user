#SJL: 10/15
#The make file for the HERCULES code using gnu compilers

#CFLAGS = -o3 -ansi -g -Wall -std=c++11 -I /Applications/C_libraries/boost_1_59_0
CFLAGS = -o3 -ansi -g -Wall -std=c++11 -I /usr/local/opt/boost/include 
LIBS = -l stdc++ -lm 
CXX = g++
OBJECTS= hercules.o initialisation.o EOS_functions.o vector_operators.o planet_functions.o concentric_layer_functions.o general_functions.o itterate.o parameters_functions.o io_routines.o

# These lines define what each file depends on, and how to make it
all: hercules #clean #add clean if needed

hercules: ${OBJECTS}
	$(CXX) $(CFLAGS) $(LIBS) ${OBJECTS} -o hercules

%.o: %.cc
	$(CXX) $(CFLAGS) $< -c -o $@
%: %.cc Makefile
	$(CXX) $(CFLAGS) $(LIBS) $< -o $@


clean:
	rm *.o
	echo Clean done

LANG = en_US.UTF-8

CPP = g++
CPPFLAGS =-std=c++11 -O3
LIBS = -lgmpxx -lgmp
OBJ = diagonalizer_q.o diagonalizer_zm.o homology_field.o matrix_zm.o clock.o
GCC_LT_4_7 := $(shell expr `g++ -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g'` \< 407)
ifeq "$(GCC_LT_4_7)" "1"
	CPPFLAGS :=-std=c++0x -O3
	LIBS += -lpthread
endif


all:$(OBJ)
	ar rcs libhomology.a diagonalizer_q.o diagonalizer_zm.o homology_field.o matrix_zm.o clock.o

%o: %.cpp 
	$(CPP) -c -o @ $< $(CPPFLAGS)

doc:
	doxygen

clean:
	rm -f $(OBJ) libhomology.a 
	rm -Rf html/ latex/

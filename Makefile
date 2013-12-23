LANG = en_US.UTF-8

CPP = g++
CPPFLAGS =-std=c++11 -O3
LIBS = -lgmpxx -lgmp
OBJ = diagonalizer_q.o diagonalizer_zm.o homology_field.o matrix_zm.o

all:$(OBJ)
	ar rcs libhomology.a diagonalizer_q.o diagonalizer_zm.o homology_field.o matrix_zm.o

%o: %.cpp 
	$(CPP) -c -o @ $< $(CPPFLAGS)

doc:
	doxygen

clean:
	rm -f $(OBJ) libhomology.a 
	rm -Rf html/ latex/

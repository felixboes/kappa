LANG = en_US.UTF-8

CPP = g++
CPPFLAGS =-std=c++11 -O3
LIBS = -lgmpxx -lgmp -lpthread
OBJ = field_coefficients.o homology_field.o clock.o

all:$(OBJ)
	ar rcs libhomology.a field_coefficients.o homology_field.o clock.o

%o: %.cpp 
	$(CPP) -c -o @ $< $(CPPFLAGS)

doc:
	doxygen

clean:
	rm -f $(OBJ) libhomology.a 
	rm -Rf html/ latex/

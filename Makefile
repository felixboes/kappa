LANG = en_US.UTF-8

CPP = g++
CPPFLAGS =-std=c++11 -O3
LIBS = -I../libhomology -lgmpxx -lgmp
OBJ = factorial.o main.o monocomplex.o tuple.o

all:$(OBJ)
	$(CPP) $(CPPFLAGS) -o monocomplex_homology $(OBJ) ../libhomology/libhomology.a $(LIBS)

%.o: %.cpp 
	$(CPP) $(CPPFLAGS) -c  $< -I../libhomology/

doc:
	doxygen

clean:
	rm -f $(OBJ) monocomplex_homology
	rm -Rf html/ latex/

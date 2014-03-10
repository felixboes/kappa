LANG = en_US.UTF-8

CPP = g++
CPPFLAGS =-std=c++11 -O3
LIBS = -I../libhomology -lgmpxx -lgmp -lboost_filesystem -lboost_system -lboost_iostreams -lboost_serialization -lpthread
OBJ = factorial.o monocomplex.o tuple.o sessionconfig.o
INCLUDES = $(wildcard *.hpp)
GCC_LT_4_7 := $(shell expr `$(CPP) -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g'` \< 407)

ifeq "$(GCC_LT_4_7)" "1"
	CPPFLAGS :=-std=c++0x -O3
endif


compute_homology:$(OBJ) ${INCLUDES} main_compute_homology.cpp
	$(CPP) $(CPPFLAGS) -o compute_homology main_compute_homology.cpp $(OBJ) ../libhomology/libhomology.a $(LIBS)

compute_statistics:$(OBJ) ${INCLUDES} main_compute_statistics.cpp
	$(CPP) $(CPPFLAGS) -o compute_statistics main_compute_statistics.cpp $(OBJ) ../libhomology/libhomology.a $(LIBS) 

compute_cache:$(OBJ) ${INCLUDES} main_compute_cache.cpp
	$(CPP) $(CPPFLAGS) -o compute_cache main_compute_cache.cpp $(OBJ) ../libhomology/libhomology.a $(LIBS) 

print_basis:$(OBJ) ${INCLUDES} main_print_basis.cpp
	$(CPP) $(CPPFLAGS) -o print_basis main_print_basis.cpp $(OBJ) ../libhomology/libhomology.a $(LIBS)

%.o: %.cpp ${INCLUDES}
	$(CPP) $(CPPFLAGS) -c  $< -I../libhomology/

doc:
	doxygen

clean:
	rm -f $(OBJ) compute_cache compute_homology compute_statistics
	rm -Rf html/ latex/

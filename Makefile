LANG = en_US.UTF-8

CPP = g++
CPPFLAGS =-std=c++11 -g
LIBS = -I../libhomology -lgmpxx -lgmp -lboost_filesystem -lboost_system -lboost_iostreams -lboost_serialization -lpthread `Magick++-config --cppflags --cxxflags --ldflags --libs`
OBJ = factorial.o monocomplex.o tuple.o sessionconfig.o css.o
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

compute_css:$(OBJ) ${INCLUDES} main_compute_css.cpp
	$(CPP) $(CPPFLAGS) -o compute_css main_compute_css.cpp $(OBJ) ../libhomology/libhomology.a $(LIBS)

draw_differentials:$(OBJ) ${INCLUDES} main_draw_differentials.cpp
	$(CPP) $(CPPFLAGS) -o draw_differentials main_draw_differentials.cpp $(OBJ) ../libhomology/libhomology.a $(LIBS)

print_basis:$(OBJ) ${INCLUDES} main_print_basis.cpp
	$(CPP) $(CPPFLAGS) -o print_basis main_print_basis.cpp $(OBJ) ../libhomology/libhomology.a $(LIBS)

%.o: %.cpp ${INCLUDES}
	$(CPP) $(CPPFLAGS) -c  $< -I../libhomology/ $(LIBS)

doc:
	doxygen

clean:
	rm -f $(OBJ) compute_cache compute_homology compute_statistics compute_css draw_differentials print_basis
	rm -Rf html/ latex/

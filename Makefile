LANG = en_US.UTF-8

CPP = g++
CPPFLAGS =-std=c++11 -O3
LIBS = -I../libhomology -lgmpxx -lgmp
OBJ = factorial.o monocomplex.o tuple.o sessionconfig.o

compute_homology:$(OBJ)
	$(CPP) $(CPPFLAGS) -o compute_homology main_compute_homology.cpp $(OBJ) ../libhomology/libhomology.a $(LIBS)

compute_statistics:$(OBJ)
	$(CPP) $(CPPFLAGS) -o compute_statistics main_compute_statistics.cpp $(OBJ) ../libhomology/libhomology.a $(LIBS) -lboost_iostreams -lboost_serialization

compute_cache:$(OBJ)
	$(CPP) $(CPPFLAGS) -o compute_cache main_compute_cache.cpp $(OBJ) ../libhomology/libhomology.a $(LIBS) -lboost_iostreams -lboost_serialization

%.o: %.cpp 
	$(CPP) $(CPPFLAGS) -c  $< -I../libhomology/

doc:
	doxygen

clean:
	rm -f $(OBJ) compute_cache compute_homology compute_statistics
	rm -Rf html/ latex/

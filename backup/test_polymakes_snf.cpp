// Tested with http://www.polymake.org/lib/exe/fetch.php/download/polymake-2.13-1.tar.bz2
//
// Put main to the polymake-2.13-1 folder and compile with
// g++ -O3 -DPOLYMAKE_DEBUG=0 -std=c++11 -o test_polymake test.cpp ./lib/core/src/Integer.cc ./lib/core/src/CharBuffer.cc ./lib/core/src/Rational.cc ./lib/core/src/inlines.cc -I./include/core/ -lgmpxx -lpthread

#include <iostream>
#include <list>

#include "polymake/Integer.h"
#include "polymake/SparseMatrix.h"
#include "polymake/Smith_normal_form.h"

int main()
{
    typedef polymake::Integer Int;
    Int a = 56;
    std::cout << a*a*a*a*a*a*a*a*a*a*a*a*a*a*a*a << std::endl;

    typedef polymake::SparseMatrix<Int> SparseMatrix;
    typedef pm::SparseMatrixStatistics<Int> SparseMatrixStatistics;
    
    SparseMatrix M(3,3);
    
    M(0,0) = 2;
    M(0,1) = 4;
    M(0,2) = 4;
    
    M(1,0) = -6;
    M(1,1) = 6;
    M(1,2) = 12;
    
    M(2,0) = 10;
    M(2,1) = -4;
    M(2,2) = -16;
    
    SparseMatrixStatistics stat;
    stat.gather(M);
    
    std::cout << stat << std::endl;
    
    std::list< std::pair<Int,int> > torsion;
    std::cout << smith_normal_form_only( M, torsion ) << std::endl;
    
    for( const auto& it : torsion )
    {
        std::cout << it.first << "; ";
    }
    std::cout << std::endl;
    
    return 0;
}


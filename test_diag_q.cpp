#include "matrix_q.hpp"
#include "diagonalizer_q.hpp"
#include <boost/numeric/ublas/io.hpp>

int main( int argc, char ** argv )
{

    MatrixQ M(4,4);
    M(0,0) = 1;
    M(0,1) = 1;
    M(0,2) = 2;
    M(0,3) = 3;

    M(1,0) = 0;
    M(1,1) = 0;
    M(1,2) = 0;
    M(1,3) = 0;

    M(2,0) = 1;
    M(2,1) = 0;
    M(2,2) = 5;
    M(2,3) = 7;

    M(3,0) = 0;
    M(3,1) = 0;
    M(3,2) = 0;
    M(3,3) = 0;

    MatrixQ N(4,2);
    N(0,0) = -15;
    N(0,1) =  -7;
    
    N(1,0) =   9;
    N(1,1) =   4;
    
    N(2,0) =   3;
    N(2,1) =   0;
    
    N(3,0) =   0;
    N(3,1) =   1;

    std::cout << "M:    " << M << std::endl;
    std::cout << "N:    " << N << std::endl;
    std::cout << "Prod: " << prod( M, N ) << std::endl;
    
    DiagonalizerQ D(M,N);
    std::cout << D.defect() << " " << D.torsion() << std::endl;
    
    return 0;
}

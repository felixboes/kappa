#include <iostream>

#include "chain_complex.hpp"
#include "diagonalizer_field.hpp"
#include "homology.hpp"
#include "matrix.hpp"
#include "matrix_zm.hpp"

int main()
{
    Zm::set_modulus(5,1);
    
    ChainComplexZm cc;
    
    ChainComplexZm::MatrixType M(4,4);
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

    ChainComplexZm::MatrixType N(4,2);
    N(0,0) = -15;
    N(0,1) =  -7;
    
    N(1,0) =   9;
    N(1,1) =   4;
    
    N(2,0) =   3;
    N(2,1) =   0;
    
    N(3,0) =   0;
    N(3,1) =   1;

    cc[2] = N;
    cc[1] = M;
    cc[0] = ChainComplexZm::MatrixType(0,4);
    
    
    cc[10] = ChainComplexZm::MatrixType(0,1);
    cc[11] = ChainComplexZm::MatrixType(1,1);
    cc(11,0,0) = 0;
    cc[12] = ChainComplexZm::MatrixType(1,1);
    cc(12,0,0) = 2;
    cc[13] = ChainComplexZm::MatrixType(0,1);
    cc[14] = ChainComplexZm::MatrixType(1,1);
    cc(14,0,0) = 2;
    cc[15] = ChainComplexZm::MatrixType(0,1);
    cc[16] = ChainComplexZm::MatrixType(1,1);
    cc(16,0,0) = 2;
    cc[17] = ChainComplexZm::MatrixType(0,1);
    cc[18] = ChainComplexZm::MatrixType(1,1);
    cc(18,0,0) = 2;
    
    
    std::cout << "N:    " << cc[2] << std::endl;
    std::cout << "M:    " << cc[1] << std::endl;
    
    ChainComplexZm::HomologyType ho = cc.homology();
    std::cout << ho << std::endl;
}

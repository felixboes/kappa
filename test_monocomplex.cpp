#include <iostream>

#include <homology.hpp>
#include "monocomplex.hpp"

template< class MatrixType >
bool equals( const MatrixType& A, const MatrixType& B )
{
    if( A.size1() == B.size1() && A.size2() == B.size2() )
    {
        uint32_t rows = A.size1();
        uint32_t cols = A.size2();
        
        for( uint32_t i = 0; i < rows; ++i )
        {
            for( uint32_t j = 0; j < cols; ++j )
            {
                if( A(i,j) != B(i,j) )
                {
                    return false;
                }
            }
        }
        
        return true;
    } 
    
    return false;
}

int main()
{
    typedef MonoComplex<ChainComplexQ> MonoComplexQ;
    MonoComplexQ mc(1,3); // p = 2( 1*2+3 )= 10
//    mc.show_basis(1);
//    mc.show_basis(2);
//    mc.show_basis(3);
//    mc.show_basis(4);
//    mc.show_basis(5);
    
//    mc.gen_differential(4);
    
    
    mc.gen_differential(1);
    mc.gen_differential_naive(1);
    mc.gen_differential(2);
    mc.gen_differential_naive(2);
    mc.gen_differential(3);
    mc.gen_differential_naive(3);
    mc.gen_differential(4);
    mc.gen_differential_naive(4);
    mc.gen_differential(5);
    mc.gen_differential_naive(5);
    mc.gen_differential(6);
    mc.gen_differential_naive(6);
    mc.gen_differential(7);
    mc.gen_differential_naive(7);
    
    
    std::cout << ( equals(mc.matrix_complex[1], mc.matrix_complex_naive[1] )? "Yeay" : "Ohh nein" ) << std::endl;
    std::cout << ( equals(mc.matrix_complex[2], mc.matrix_complex_naive[2] )? "Yeay" : "Ohh nein" ) << std::endl;
    std::cout << ( equals(mc.matrix_complex[3], mc.matrix_complex_naive[3] )? "Yeay" : "Ohh nein" ) << std::endl;
    std::cout << ( equals(mc.matrix_complex[4], mc.matrix_complex_naive[4] )? "Yeay" : "Ohh nein" ) << std::endl;
    std::cout << ( equals(mc.matrix_complex[5], mc.matrix_complex_naive[5] )? "Yeay" : "Ohh nein" ) << std::endl;
    std::cout << ( equals(mc.matrix_complex[6], mc.matrix_complex_naive[6] )? "Yeay" : "Ohh nein" ) << std::endl;
    std::cout << ( equals(mc.matrix_complex[7], mc.matrix_complex_naive[7] )? "Yeay" : "Ohh nein" ) << std::endl;
    
    return 0;
}

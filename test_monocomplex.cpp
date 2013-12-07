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

void print_usage(int argc, char** argv)
{
    std::cout << "Usage: " << argv[0] << " g m differential_to_print (show_basis = 0 or 1) (check_dd_zero = 0 or 1)" << std::endl;
}

int main(int argc, char** argv)
{
    if( argc < 4 )
    {
        print_usage(argc, argv);
        return 1;
    }
    
    typedef MonoComplex<ChainComplexQ> MonoComplexQ;
    MonoComplexQ mc( atoi(argv[1]), atoi(argv[2]) ); // p = 2( 1*2+3 )= 10
    
    if( argc >= 5 && atoi(argv[4]) != 0 )
    {
        mc.show_basis( atoi(argv[3]) - 1 );
        mc.show_basis( atoi(argv[3]) );
    }
    
    mc.gen_differential( atoi(argv[3]) );
    
    std::cout << mc.matrix_complex[atoi(argv[3])] << std::endl;
    
    if( argc >= 6 && atoi(argv[5]) != 0 )
    {
        mc.gen_differential( atoi(argv[3]) - 1 );
        std::cout << mc.matrix_complex[atoi(argv[3]) - 1] << std::endl;
        
        std::cout << "Check if dd = 0:" << std::endl;
        std::cout << prod( mc.matrix_complex[atoi(argv[3]) - 1], mc.matrix_complex[atoi(argv[3])] ) << std::endl; 
    }
    
    return 0;
}

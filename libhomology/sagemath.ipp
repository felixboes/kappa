#include "sagemath.hpp"

template< class CoefficientType >
void SagemathInterface::create_matrix( size_t rows, size_t cols )
{
    std::stringstream command;
    command << "del A\nA=matrix( " << sage_coeff<CoefficientType>() << ", " << rows << ", " << cols << ", sparse=True)\n";
    fprintf( sagemath_pipe, command.str().c_str() );
}

template< class MatrixType >
void SagemathInterface :: update_row( const int32_t row, const MatrixType& vect )
{
    if( sagemath_pipe == nullptr )
    {
        std::cout << "pipe closed" << std::endl;
        return;
    }

    if( vect.size2() == 0 )
    {
        return;
    }

    std::stringstream command;
    command << "A.set_row(" << row << ", [" << vect.at(0,0);
    for( size_t j = 1; j < vect.size2(); ++j )
    {
        command << "," << vect.at(0,j);
    }
    command << "])\n";

    fprintf( sagemath_pipe, command.str().c_str() );
}

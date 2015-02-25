#include "kappa.hpp"

void play_in_the_double_complex()
{
    Tuple::parallel_case();
    
    DoubleComplexBasis basis;
    
    basis.add_basis_element( create_cell(2, 4, 3, 2, 1) );
    basis.add_basis_element( create_cell(2, 4, 1, 3, 2) );
    basis.add_basis_element( create_cell(2, 2, 1, 4, 3) );
    basis.add_basis_element( create_cell(2, 3, 2, 4, 1) );
    
    for( const auto& cell : basis.basis )
    {
        std::cout << cell << std::endl;
        
        const std::map< uint8_t, int8_t > or_sign( cell.orientation_sign() );
    
        for( uint32_t i = Tuple::get_min_boundary_offset(); i <= 4 - Tuple::get_max_boundary_offset(); i++ )
        {
            Tuple boundary;
            if( (boundary = cell.d_hor_double_complex(i)) )
            {
                std::cout << "    " << ( i +or_sign.at(i) % 2 == 0 ? " 1 " : "-1 " ) << boundary << std::endl;
            }
        }
    }
}

int main( int argc, char** argv )
{
    std::cout.setf(std::ios::unitbuf);
    
    std::cout << kappa_version( argc, argv ) << std::endl;
    
    if( argc > 1 )
    {
        if( atoi( argv[1] ) == 1 )
        {
            play_in_the_double_complex();
            return 0;
        }
    }
    
    play_in_the_double_complex();
    return 0;
    
}

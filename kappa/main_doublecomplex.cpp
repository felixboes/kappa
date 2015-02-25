#include <boost/range/adaptor/reversed.hpp>

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
                std::cout << "    " << ( (1-2*(i%2)) * or_sign.at(i) == 1 ? " 1 " : "-1 " ) << boundary << std::endl;
            }
        }
    }
}


void create_and_test_doublecomplex(uint32_t g = 0, uint32_t m = 2)
{
    std::cout << "Rational: " << std::endl;
    DoubleComplex< ChainComplexQ > dcq(g, m, SignConvention::all_signs, 1, 1);
    
    for( const auto& it : boost::adaptors::reverse( dcq.bases ) )
    {
        std::cout << "Topcells of horizontal dimension " << (int32_t)it.first << ":" << std::endl;
        std::cout << it.second << std::endl;
    }
    
    dcq.gen_differential(4);
    
    std::cout << dcq.get_current_differential() << std::endl;
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
        else if( atoi( argv[1] ) == 2 )
        {
            create_and_test_doublecomplex();
            return 0;
        }
    }
    
    create_and_test_doublecomplex();
    
    return 0;
    
}

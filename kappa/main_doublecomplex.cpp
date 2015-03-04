#include <boost/range/adaptor/reversed.hpp>

#include "kappa.hpp"

void play_in_the_double_complex()
{   
    DoubleComplexBasis basis;
    
    basis.add_basis_element( create_highcell(2, 4, 3, 2, 1) );
    basis.add_basis_element( create_highcell(2, 4, 1, 3, 2) );
    basis.add_basis_element( create_highcell(2, 2, 1, 4, 3) );
    basis.add_basis_element( create_highcell(2, 3, 2, 4, 1) );
    
    for( const auto& cell : basis.basis_col )
    {
        std::cout << cell << std::endl;
        
        const std::map< uint8_t, int8_t > or_sign( cell.orientation_sign() );
    
        for( uint32_t i = HighCell::get_min_boundary_offset(); i <= 4 - HighCell::get_max_boundary_offset(); i++ )
        {
            HighCell boundary;
            if( (boundary = cell.d_hor_double_complex(i)) )
            {
                std::cout << "    " << ( (1-2*(i%2)) * or_sign.at(i) == 1 ? " 1 " : "-1 " ) << boundary << std::endl;
            }
        }
    }
    
    for( const auto& cell : basis.basis_ess )
    {
        std::cout << cell << std::endl;
        
        const std::map< uint8_t, int8_t > or_sign( cell.orientation_sign() );
    
        for( uint32_t i = HighCell::get_min_boundary_offset(); i <= 4 - HighCell::get_max_boundary_offset(); i++ )
        {
            HighCell boundary;
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
        std::cout << "Cells of horizontal degree " << (int32_t)it.first << ":" << std::endl;
        std::cout << it.second << std::endl;
    }
    
    const size_t h = 2*g+m;
    for( size_t p = 2*h; p >= h; --p )
    {
        dcq.gen_differential(p);
        auto& diff = dcq.get_current_differential();
        std::cout << diff << std::endl;
        DiagonalizerField<MatrixQ> diago;
        diago.diag_field(diff);
        std::cout << diff << std::endl;
        
        dcq.compute_proj_E(p);
        std::cout << diff << std::endl;
    }
}

int main( int argc, char** argv )
{
    std::cout.setf(std::ios::unitbuf);
    
    std::cout << kappa_version( argc, argv ) << std::endl;
    
    if( argc > 2 )
    {
        create_and_test_doublecomplex( atoi(argv[1]) , atoi( argv[2] ) );
        return 0;
    }
    
    create_and_test_doublecomplex();
    
    return 0;
    
}

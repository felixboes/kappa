#include <boost/range/adaptor/reversed.hpp>

#include "kappa.hpp"

void play_in_the_double_complex()
{   
    DoubleComplexBasis basis;
    
    basis.add_basis_element( create_highcell(2, 4, 3, 2, 1) );
    basis.add_basis_element( create_highcell(2, 4, 1, 3, 2) );
    basis.add_basis_element( create_highcell(2, 2, 1, 4, 3) );
    basis.add_basis_element( create_highcell(2, 3, 2, 4, 1) );
    
    basis.generate_indices();
    
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
    
    std::cout << basis << std::endl;
    
    const Tuple t = create_cell(2,4,3,2,1);
    std::cout << t << " " << basis.id_of(t) << std::endl;
}

void create_and_test_doublecomplex(uint32_t g = 1, uint32_t m = 1)
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

void test_proj_E()
{
    DoubleComplex< ChainComplexQ > dcq(0, 2, SignConvention::all_signs, 1, 1);
    dcq.compute_proj_E(4);
    dcq.proj_E_ast( create_cell(2, 3, 2, 2, 1) );
    dcq.proj_E_ast( create_cell(2, 3, 1, 3, 2) );
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
    
    test_proj_E();
    
    return 0;
    
}

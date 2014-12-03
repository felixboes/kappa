#include <iostream>

#include "homology.hpp"
#include <kappa/kappa.hpp>

int main(int argc, char** argv)
{
    if( argc == 1 )
    {
        std::cout << "GMP Version: " << gmp_version << std::endl;
        
        VectorQ v(3);
        v(0) = 1;
        v(1) = 0;
        v(2) = 0;
        
        // 2,  1,  1
        // 1, -1,  2
        // 0,  1, -1
        MatrixQ m(3,3);
        m(0,0) = 2;
        m(0,1) = 1;
        m(0,2) = 1;
        m(1,0) = 1;
        m(1,1) = -1;
        m(1,2) = 2;
        m(2,0) = 0;
        m(2,1) = 1;
        m(2,2) = -1;
        
        std::cout << m << std::endl;
        DiagonalizerField<MatrixQ> diago;
        diago(m);
        std::cout << m << std::endl;
        std::cout << m.base_changes() << std::endl;;
        std::cout << m.triangular_shape() << std::endl;
        
        VectorQ v_prime(3);
        v_prime(0) =  1;
        v_prime(1) = -1;
        v_prime(2) = -1;
        apply_base_changes_kernel( m, v_prime );
        std::cout << v_prime << std::endl;
        
        v_prime(0) = 0;
        v_prime(1) = 1;
        v_prime(2) = 0;
        apply_base_changes_image( m, v_prime );
        std::cout << v_prime << std::endl;
        
        
        VectorQ v_2(3);
        v_2(0) = 3;
        v_2(1) = Q(7) / Q(4);
        v_2(2) = -2;
        
        v += v_2;
        std::cout << v << std::endl;
    
        m.cache_matrix("test_matrix");
        std::cout << load_from_file_bz2< MatrixField<Q> > ("test_matrix") << std::endl;
        
        save_to_file_bz2( m.base_changes(), "test_matrix" );
        std::cout << load_from_file_bz2< MatrixField<Q> > ("test_matrix") << std::endl;
        
        save_to_file_bz2( m.triangular_shape(), "test_matrix" );
        std::cout << load_from_file_bz2< MatrixField<Q> > ("test_matrix") << std::endl;
        /*
        // 0,  1,  1
        // 1, -1,  0
        // 0,  1, -1
        MatrixBool n(3,3);
        n.add_entry(0,1);
        n.add_entry(0,2);
        n.add_entry(1,0);
        n.add_entry(1,1);
        n.add_entry(2,1);
        n.add_entry(2,2);
        
        VectorBool w(3);
        w.add_entry(0);
        w.add_entry(2);
        
        std::cout << n << std::endl;
        DiagonalizerField<MatrixBool> diagobool;
        diagobool(n);
        std::cout << n << std::endl;
        n.print_triangular_shape();
        n.print_base_changes_in_short_form();
        
        VectorBool w_prime(3);
        w_prime.add_entry(1);
        w_prime.add_entry(2);
        
        apply_base_changes<MatrixBool, VectorBool>( n, w_prime );
        std::cout << w_prime << std::endl;
        
        VectorBool w_2(3);
        w_2.add_entry(0);
        w_2.add_entry(1);
        
        std::cout << w << " + " << w_2 << " = ";
        w += w_2;
        std::cout << w << std::endl;
        
        v.resize(10);
        std::cout << v << std::endl;
        
        w.resize(10);
        std::cout << w << std::endl;
        */
    }
    else
    {
        std::cout << load_from_file_bz2< MatrixField<Q> >( argv[1] );
    }
    return 0;
}

#include "generators.hpp"

template< class CoefficientT >
MonoCochainField< CoefficientT > create_cochain( const Generator& name )
{
    typedef MonoCochainField< CoefficientT > CochainType;
    switch( name )
    {
    
    case a:
    {
        CochainType cochain(0, 1, 2);
        cochain.set_name("a");
        cochain( create_cell(1, 2, 1) ) = CoefficientT(1);
        return cochain;
    }
    case b:
    {
        CochainType cochain(0, 2, 3);
        cochain.set_name("b");
        cochain( create_cell(2, 3, 2, 2, 1) ) = CoefficientT(1);
        cochain( create_cell(2, 3, 1, 3, 2) ) = CoefficientT(1);
        return cochain;
    }
    case c:
    {
        CochainType cochain(1, 0, 4);
        cochain.set_name("c");
        cochain( create_cell(2, 4, 2, 3, 1 ) ) = CoefficientT(1);
        return cochain;
    }
    case d:
    {
        CochainType cochain(1, 0, 3);
        cochain.set_name("d");
        cochain( create_cell(2, 3, 1, 2, 1 ) ) = CoefficientT(1);
        return cochain;
    }
    case e:
    {
        CochainType cochain(1, 1, 4);
        cochain.set_name("e");
        cochain.add_kappa_dual( CoefficientT(1), create_cell(3, 3, 1, 4, 3, 2, 1) );
        return cochain;
    }
    case Qc:
    {
        CochainType cochain(2, 0, 7);
        cochain.set_name("Qc");
        // Observe that every two adjecent cells differ by a face index 1, hence the sign is kept.
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 7, 5, 6, 4, 4, 2, 3, 1) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 7, 5, 6, 1, 4, 2, 3, 1) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 7, 5, 6, 2, 4, 2, 3, 1) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 7, 5, 6, 3, 4, 2, 3, 1) );
        // Now jumping to the groud differs by a face index of 4, hence we have to change the sign.
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 5, 6, 1, 5, 3, 4, 2) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 2, 6, 1, 5, 3, 4, 2) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 3, 6, 1, 5, 3, 4, 2) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 4, 6, 1, 5, 3, 4, 2) );
        // Now jumping to the groud differs by a face index of 4, hence we have to change the sign.
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 7, 2, 6, 1, 6, 4, 5, 3) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 7, 2, 3, 1, 6, 4, 5, 3) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 7, 2, 4, 1, 6, 4, 5, 3) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 7, 2, 5, 1, 6, 4, 5, 3) );
        // Now jumping to the groud differs by a face index of 4, hence we have to change the sign.
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 2, 3, 1, 7, 5, 6, 4) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 4, 2, 3, 1, 7, 5, 6, 4) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 5, 2, 3, 1, 7, 5, 6, 4) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 6, 2, 3, 1, 7, 5, 6, 4) );
        
        return cochain;
    }
    case Qd:
    {
        CochainType cochain(2, 0, 5);
        cochain.set_name("Qd");
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 3, 3, 1, 2, 1 ) );
    
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 1, 3, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 4, 1, 3, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 2, 3, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 4, 2, 3, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 2, 4, 2, 3, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 4, 1, 4, 2, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 2, 4, 1, 4, 2, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 1, 4, 2, 3, 2 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 2, 1, 4, 2, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 2, 1, 4, 2, 3, 2 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 3, 1, 4, 2, 3, 2 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 2, 1, 5, 3, 4, 3 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 3, 1, 2, 1, 5, 3, 4, 3 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 4, 1, 2, 1, 5, 3, 4, 3 ) );
        
        return cochain;
    }
    case Te:
    {
        CochainType cochain(2, 0, 5);
        cochain.set_name("Te");
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 5, 3, 3, 1, 4, 3, 2, 1) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 5, 1, 3, 1, 4, 3, 2, 1) );
        return cochain;
    }
        
    }
    
    
    return CochainType(0,1,2);
}



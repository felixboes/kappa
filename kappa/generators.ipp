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
        
    }
    
    
    return CochainType(0,1,2);
}



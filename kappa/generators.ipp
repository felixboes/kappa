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
    }
    
    return CochainType(0,1,2);
}



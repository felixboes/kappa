#include "homology_field.hpp"

HomologyField::HomologyField()
{
}

HomologyField::HomologyField( int32_t n, KernT ker, TorsT tor )
{
    kern[n] = ker;
    tors[n] = tor;
}

void HomologyField::set_kern( int32_t n, KernT k )
{
    kern[n] = k;
}

void HomologyField::set_tors( int32_t n, TorsT t )
{
    tors[n] = t;
}

std::ostream &operator << ( std::ostream& os, const HomologyField& homol )
{
    for( auto it = homol.kern.begin(); it != homol.kern.end(); ++it )
    {
        int64_t torsion;
        if( homol.tors.count( it->first ) == 0 )
        {
            torsion = 0;
        }
        else
        {
            torsion = homol.tors.at( it->first );
        }
        
        os << "Homologiemodul an der Stelle " << std::setw(4) << it->first << std::endl;
        os << std::setfill('-') << std::setw(40) << "-" << std::setfill(' ') << std::endl;
        os << "Dimension = " << it->second - torsion << std::endl;
        os << std::endl;
    }
    return os;
}

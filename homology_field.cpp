#include "homology_field.hpp"

HomologyField::HomologyField()
{
}

HomologyField::HomologyField( int32_t n, KernT kern, TorsT tors )
{
    dimension[n] = kern - tors;
}

std::ostream &operator << ( std::ostream& os, const HomologyField& homol )
{
    for( auto& it : homol.dimension )
    {
        os << "Homologiemodul an der Stelle " << std::setw(4) << it.first << std::endl;
        os << std::setfill('-') << std::setw(40) << "-" << std::setfill(' ') << std::endl;
        os << "Dimension = " << it.second << std::endl;
        os << std::endl;
    }
    return os;
}

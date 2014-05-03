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

HomologyField::KernT HomologyField::get_kern( int32_t n ) const
{
    if( kern.count(n) == 0 )
    {
        return KernT();;
    }
    else
    {
        return kern.at(n);
    }
}

HomologyField::TorsT HomologyField::get_tors( int32_t n ) const
{
    if( tors.count(n) == 0 )
    {
        return TorsT();
    }
    else
    {
        return tors.at(n);
    }
}

void HomologyField::erase_kern( int32_t n )
{
    kern.erase(n);
}

void HomologyField::erase_tors( int32_t n )
{
    tors.erase(n);
}

std::ostream &operator << ( std::ostream& os, const HomologyField& homol )
{
    for( const auto& it : homol.kern )
    {
        auto & p = it.first;
        auto & kern = it.second;
        
        os << "Homology module H_" << p << std::endl;
        os << std::setfill('-') << std::setw(35) << "-" << std::setfill(' ') << std::endl;
        os << "Dimension = " << kern - homol.get_tors(p) << std::endl;
        os << std::endl;
    }
    return os;
}

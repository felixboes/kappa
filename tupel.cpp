#include "tupel.hpp"

Tupel :: Tupel(size_t h) :
    p(0),
    rep( h, std::pair< uint8_t, uint8_t >(0, 0) )
{
}

Transposition& Tupel :: operator[](size_t n)
{
    return rep[n];
}

int32_t Tupel :: norm() const
{
    return rep.size();
}

bool Tupel :: valid() const
{
    for( auto const &it : rep )
    {
        if( it.first == 0 || it.second == 0 )
        {
            return false;
        }
    }
    
    return true;
}

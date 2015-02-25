#include "monobasis.hpp"

MonoBasis::MonoBasis() : basis()
{
}
 
uint32_t MonoBasis :: add_basis_element (Tuple t)
{
    uint32_t id = basis.size();
    t.id = id;
    basis.insert(std::move(t));
    
    return id;
}

uint64_t MonoBasis :: size() const
{
    return basis.size();
}

int64_t MonoBasis :: id_of(const Tuple &t) const
{
    auto it = basis.find(t);
    if( it == basis.end() )
    {
        return -1;
    }
    else
    {
        return it->id;
    }
}

MonoBasis load_parallel_mono_basis( const uint32_t g, const uint32_t m, const int32_t p )
{
    std::string filename =
            "./cache/bases_parallel/" +
            std::to_string(g) + "_" +
            std::to_string(m) + "_" +
            std::to_string(p);
    
    return load_from_file_bz2< MonoBasis >( filename, false );
}

std::ostream& operator<< ( std::ostream& os, const MonoBasis& mb )
{
    for( const auto it : mb.basis )
    {
        os << it.id << ": " << it << std::endl;
    }
    return os;
}

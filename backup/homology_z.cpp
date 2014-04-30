#include "homology_z.hpp"

homology::homology()
{
}

uint32_t homology::free_part()
{
    return free_part;
}

uint32_t tors_part()
{
    return tors_part;
}

Coefficient &operator[] (uint32_t n)
{
    return tors_[num]
}

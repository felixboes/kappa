#ifndef HOMOLOGY_DUMMY_HPP
#define HOMOLOGY_DUMMY_HPP

#include <iostream>

/**
 *  Dummy class. See HomologyField
 */
class HomologyDummy
{
public:
    typedef int64_t KernT;
    typedef int64_t TorsT;
    HomologyDummy() {}
    HomologyDummy( int32_t, KernT, TorsT ) {}
    inline void set_kern( int32_t n, KernT k ) {}
    inline void set_tors( int32_t n, TorsT t ) {}
    friend std::ostream& operator<< (std::ostream& stream, const HomologyDummy& homol);
};

inline std::ostream& operator<< (std::ostream& stream, const HomologyDummy& homol)
{
    return stream << "This is a dummy module for homology computations.";
}

#endif // HOMOLOGY_DUMMY_HPP

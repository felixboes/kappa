#ifndef HOMOLOGY_DUMMY_HPP
#define HOMOLOGY_DUMMY_HPP

// Description:
//
// This header defines a dummy version of homology.
// This is used if we are not at all interested in homology computations, for example when we construct and save differentials.

#include <iostream>

/**
 *  This is a dummy class that mimes the functionality of HomologyField, but does nothing.
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

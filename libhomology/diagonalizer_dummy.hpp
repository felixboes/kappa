#ifndef DIAGONALIZER_DUMMY_HPP
#define DIAGONALIZER_DUMMY_HPP

// Description:
//
// This header defines a dummy version of a diagonalizer.
// This is used if we are not at all interested in homology computations, for example when we construct and save differentials.

#include <boost/numeric/ublas/io.hpp>

#include "homology_dummy.hpp"
#include "parallelization.hpp"

/**
 *  This is a dummy class that mimes the functionality of DiagonalozerField, but does nothing.
 */
template <class MatrixClass>
class DiagonalizerDummy
{
public:
    DiagonalizerDummy() {}
    void operator() ( MatrixClass &matrix, uint32_t number_threads=0 ) {}
    void operator() ( MatrixClass &matrix, atomic_uint & current_rank, uint32_t number_threads=0 ) {}
    uint32_t dfct() {return 0;}
    HomologyDummy::KernT kern() {return HomologyDummy::KernT();}
    uint32_t rank() {return 0;}
    HomologyDummy::TorsT tors() {return HomologyDummy::TorsT();}
};

#endif // DIAGONALIZER_DUMMY_HPP

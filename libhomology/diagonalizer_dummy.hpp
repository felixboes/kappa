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
template <class MatrixType>
class DiagonalizerDummy
{
public:
    DiagonalizerDummy() : transp(false) {}
    void operator() ( MatrixType & ) {}
    uint32_t dfct() {return 0;}
    HomologyDummy::KernT kern() {return HomologyDummy::KernT();}
    uint32_t rank() {return 0;}
    HomologyDummy::TorsT tors() {return HomologyDummy::TorsT();}
    
    void apply_base_changes( MatrixType& differential, const MatrixType& base_changes ){ (void)differential; (void)base_changes; } // Supress 'unused variable' warnings.

    bool transp;
    uint32_t num_working_threads;
    uint32_t num_remaining_threads;
};

#endif // DIAGONALIZER_DUMMY_HPP

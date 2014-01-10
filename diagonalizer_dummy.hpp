#ifndef DIAGONALIZER_DUMMY_HPP
#define DIAGONALIZER_DUMMY_HPP

#include <boost/numeric/ublas/io.hpp>

#include "homology_dummy.hpp"
#include "parallelization.hpp"

/**
 *  Dummy class. Compare with DiagonalizerZm.
 */

template <class MatrixClass>
class DiagonalizerDummy
{
public:
    DiagonalizerDummy() {}
    void operator() ( MatrixClass &matrix, uint32_t number_threads=0 ) {}
    void operator() ( MatrixClass &matrix, atomic_uint & current_rank, uint32_t number_threads=0 ) {}
    void operator() ( MatrixClass &post_matrix, MatrixClass &matrix ) {}
    void operator() ( MatrixClass &post_matrix, MatrixClass &matrix, MatrixClass &pre_matrix );
    uint32_t dfct() {return 0;}
    HomologyDummy::KernT kern() {return HomologyDummy::KernT();}
    uint32_t rank() {return 0;}
    HomologyDummy::TorsT tors() {return HomologyDummy::TorsT();}
};

#endif // DIAGONALIZER_DUMMY_HPP

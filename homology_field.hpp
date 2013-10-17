#ifndef HOMOLOGY_FIELD_HPP
#define HOMOLOGY_FIELD_HPP

#include <cinttypes>
#include <iomanip>
#include <iostream>
#include <map>

/**
 * Class implementing the homology of a chain complex with coefficients in a field.
 * Over a field, a homology module is characterized by the defect of the
 * outgoing differential and the dimension of the rank of the ingoing differential.
 */
class HomologyField
{
public:
    /// type of the kernel
    typedef int64_t KernT;

    /// type of the torsion
    typedef int64_t TorsT;
    
    /// trivial constructor
    HomologyField();

    /// constructor initializing a single homology module
    HomologyField( int32_t, KernT, TorsT );
    
    /// sets the kernel
    void set_kern( int32_t n, KernT k );

    /// sets the torsion
    void set_tors( int32_t n, TorsT t );
    
    /// output stream
    friend std::ostream& operator<< (std::ostream& stream, const HomologyField& homol);

private:
    /// stores the kernel for each possibly non-zero homology module
    std::map< int32_t, int64_t > kern;

    /// stores the torsion for each possibly non-zero homology module
    std::map< int32_t, int64_t > tors;
};

std::ostream& operator<< (std::ostream& stream, const HomologyField& homol);

#endif // HOMOLOGY_FIELD_HPP

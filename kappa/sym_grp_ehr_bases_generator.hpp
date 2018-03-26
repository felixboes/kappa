#ifndef KAPPA_SYM_GRP_EHR_BASES_GENERATOR_H
#define KAPPA_SYM_GRP_EHR_BASES_GENERATOR_H

#include "sym_grp_tuple.hpp"
#include "ehr_basis.hpp"

/**
 * This class generates the bases of the modules in the Ehrenfried complex wrt.:
 *     genus g,
 *     number of punctures m,
 *     radial resp. parallel case according to SymGrpTuple::is_radial(),
 *     the symmetric group, ie. all entries in the tuples are Transpositions.
 */
class SymGrpEhrBasesGenerator
{

public:

    SymGrpEhrBasesGenerator(uint8_t _g, uint8_t _m);

    std::map< int32_t, EhrBasis >& generate_bases();


private:

    /** Recursive function initializing the basis_complex.
        In the call of gen_bases with the parameters s, p and tuple, we assume that the first s transpositions
        containing symbols 1, ..., p are fixed and append all possible transpositions at position s+1, applying
        the function recursively in an appropriate way.
        If s == h, we don't append another transposition since we have completed a possible basis element.
        We check whether its number of cycles is appropriate for it to be a basis element, and if this is
        the case, we add it to the basis in degree p. Thereby, basis elements are sorted according to the number
        of clusters.
    **/
    void gen_bases( const uint32_t s, const uint32_t p, const uint32_t start_symbol, SymGrpTuple& tuple);

    uint8_t g;
    uint8_t m;
    uint8_t h;

    // basis_complex[p] is the basis of the p-th module in the EhrComplex associated to g, m, radial.
    // computed only when generate_bases() is called beforehand.
    std::map< int32_t, EhrBasis > basis_complex;
};


#endif //KAPPA_SYM_GRP_EHR_BASES_GENERATOR_H

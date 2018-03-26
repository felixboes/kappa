#ifndef KAPPA_ALT_GRP_EHR_BASES_GENERATOR_H
#define KAPPA_ALT_GRP_EHR_BASES_GENERATOR_H

#include "ehr_basis.hpp"
#include "alt_grp_tuple.hpp"

/**
 * This class generates the bases of the modules in the Ehrenfried complex wrt.:
 *     genus g,
 *     number of punctures m,
 *     radial resp. parallel case according to AltGrpTuple::is_radial(),
 *     the alternating group, ie. all entries in the tuples are Norm2Permutations.
 */
class AltGrpEhrBasesGenerator
{

public:

    AltGrpEhrBasesGenerator(uint8_t _g, uint8_t _m);

    // only possible if h is even, where h = 2*_g+_m-1 + AlgGrpTuple::get_min_symbol().
    std::vector< std::vector< AltGrpTuple> >& generate_bases();

    void print_bases();



// PUBLIC ONLY FOR UNIT TESTS

    void generate_possible_norm2permutations();

    /**
     * @param curr_tuple: tuple with h/2 entries, where the first k:=curr_num_entries are already set s.t. curr_tuple is
     * unstable in these first k entries and s.t. the maximum symbol used in the first k entries is curr_max_symbols.
     *
     * @note: k=0 is allowed.
     *
     * If k<h/2, computes all possibilities for the (k+1)-th entry st. the tuple is unstable at this position (k, k+1).
     * If k=h/2, adds curr_tuple to the according basis if it is a cell.
     */
    void generate_bases_recursively(uint8_t curr_num_entries, uint8_t curr_max_symbol, AltGrpTuple& curr_tuple);


    uint8_t g;
    uint8_t m;
    uint8_t h;

    // possible_norm2permutations[i] are all norm2permutations with symbols in {min_symbol,...,2h} with height i<=2h,
    // written in norm2notation
    std::vector< std::vector< Norm2Permutation> > possible_norm2permutations;

    // bases[p] is the basis of the p-th module in the EhrenfriedComplex associated to g, m, radial.
    // computed only when generate_bases() is called beforehand.
    std::vector< std::vector< AltGrpTuple > > bases;
};


#endif //KAPPA_ALT_GRP_EHR_BASES_GENERATOR_H

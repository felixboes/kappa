#include "sym_grp_ehr_bases_generator.hpp"

SymGrpEhrBasesGenerator::SymGrpEhrBasesGenerator(uint8_t _g, uint8_t _m) :
    g(_g),
    m(_m)
{
    h=2*g+m;
    if(SymGrpTuple::is_radial())
    {
        h--;
    }
}

void SymGrpEhrBasesGenerator::gen_bases(const uint32_t s, const uint32_t p, const uint32_t start_symbol,
                                        SymGrpTuple &tuple)
{
    /* Up to now we have determined all fully unstable tuples of s transpositions containing the
       symbols 1, ..., p, each at least once. We now add an (s+1)-th transposition and continue
       recursively.*/
    if(s < h) // There are h-s transpositions left to be determined.
    {
        /* From an s-tuple containing p symbols we can build up an (s+1)-tuple
        with p, p+1 or p+2 symbols. */

        /* p -> p
           In this case we use the same number of symbols. Since we only
           enumerate fully unstable tuples, the height of the (s+1)-th transposition needs to be p.
           We try out all possibilities for the second symbol in the (s+1)-th transposition. */
        tuple.p = p;
        for(uint32_t i = start_symbol; i < p; ++i)
        {
            tuple[s+1] = Transposition(p, i);
            gen_bases(s+1, p, start_symbol, tuple);
        }

        /* p -> p+1
           In this case we use p+1 symbols instead of p. Two cases occur. */
        tuple.p = p+1;
        /* Case 1: The new row in the parallel slit domain is inserted at the top, i.e.
                   the transpositions 1, ..., s remain the same and we only insert the symbol
                   p+1 in the (s+1)-th transposition, together with any symbol of 1, ..., p. */
        for (uint32_t i = start_symbol; i <= p; ++i)
        {
            tuple[s+1] = Transposition(p+1, i);
            gen_bases(s+1, p+1, start_symbol, tuple);
        }

        /* Case 2: The new row in the parallel slit domain is not inserted at the top but at a position
                   i = 1, ..., p. Thus all indices i+1, ..., p are shifted up by one. */
        for(uint32_t i = start_symbol; i <= p; ++i)
        {
            SymGrpTuple tmp = tuple;
            for( uint32_t j = s; j >= 1; --j )
            {
                if( tmp[j].first >= i )
                {
                    tmp[j].first++;
                    if(tmp[j].second >= i )
                    {
                        tmp[j].second++;
                    }
                }
                else
                {
                    break;
                }
            }

            tmp[s+1] = Transposition(p+1, i);
            gen_bases(s+1, p+1, start_symbol, tmp);
        }

        /* p -> p+2
           Now we use p+2 symbols instead of p. Thus one row is inserted at the top of the
           parallel slit domain, i.e. p+2 is the height of the (s+1)-th transposition. The other row is inserted
           either directly below the top row or between the old rows. Since the symbol p+1 does not
           occur in the transpositions 1, ..., s, both cases can be expressed by choosing a symbol
           i = 1, ..., p+1 and by shifting up all indices >= i by one. */
        tuple.p = p+2;
        for( uint32_t i = start_symbol; i <= p + 1; ++i)
        {
            SymGrpTuple tmp = tuple;
            for( uint32_t j = s; j >= 1; j-- )
            {
                if( tmp[j].first >= i )
                {
                    tmp[j].first++;
                    if(tmp[j].second >= i )
                    {
                        tmp[j].second++;
                    }
                }
                else
                {
                    break;
                }
            }

            tmp[s+1] = Transposition(p+2, i);
            gen_bases(s+1, p+2, start_symbol, tmp);
        }
    }
    else // Check whether the created h-tuple is really a generator, i.e. if it has the correct
         // number of cycles. If this is the case, we add tuple to the basis elements of the
         // p-th basis and store the index of tuple in this basis as the id of tuple.
    {
        if (tuple.has_correct_num_cycles(m))
        {
            tuple.id = basis_complex[p].add_basis_element( tuple );
        }
    }
}

std::map< int32_t, EhrBasis >& SymGrpEhrBasesGenerator::generate_bases()
{
    // Generate all tuples with h transpositions containing the symbols 1, ..., p,
    // each at least once, with the correct number of cycles.
    SymGrpTuple tuple(h);
    tuple[1] = Transposition(2, 1);
    tuple.p = 2;
    gen_bases(1, 2, 1, tuple);  // We start with the transposition ... (2 1).
    // In the radial case, we also generate all tuples as above, but additionally all tuples containing
    // the symbol 0.
    if (SymGrpTuple::is_radial() )
    {
        SymGrpTuple radial_tuple(h);
        radial_tuple.p = 1;
        radial_tuple[1] = Transposition(1, 0);
        gen_bases(1, 1, 0, radial_tuple);
    }
    return basis_complex;
}

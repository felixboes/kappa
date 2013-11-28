#include "monocomplex.hpp"

/*
 *
 *   MonoComplex
 *
 */
template< class MatrixComplex >
MonoComplex< MatrixComplex > :: MonoComplex(uint32_t _g, uint32_t _m) : g(_g), m(_m), h(2*_g + _m)
{
    Tuple tuple(h);
    tuple[1] = Transposition(2, 1);
    gen_bases(1, 2, tuple);  // We start with the transposition ... (2 1).
}

template< class MatrixComplex >
void MonoComplex< MatrixComplex > :: show_basis( int32_t p ) const
{
    if( basis_complex.count(p) )
    {
         std::cout << "This it the " << p << "-th basis:" << std::endl;
        const auto& basis_vector = basis_complex.at(p).basis;
        for( auto it = basis_vector.cbegin(); it != basis_vector.cend(); ++it )
        {
            std::cout << it->id << ": " << *it << std::endl;
        }
    }
    else
    {
        std::cout << "The " << p << "-th basis is empty" << std::endl;
    }
}

template< class MatrixComplex >
void MonoComplex< MatrixComplex > :: gen_bases(uint32_t l, uint32_t p, Tuple& tuple)
{
    /* Up to now we have determined all monotonic tuples of l transpositions containing the 
       symbols 1, ..., p, each at least once. We now add an (l+1)-th transposition and continue
       recursively.*/
    if(l < h) // There are h-l transpositions left to be determined.
    {
	/* From an l-tuple containing p symbols we can build up an (l+1)-tuple 
        with p, p+1 or p+2 symbols. */

        /* p -> p
           In this case we use the same number of symbols. Since we only
           enumerate monotonic tuples, the height of the (l+1)-th transposition needs to be p.
           We try out all possibilities for the second symbol in the (l+1)-th transposition. */
        tuple.p = p;
        for(uint32_t i = p-1; i > 0; i--)
        {
            tuple[l+1] = Transposition(p, i); 
            gen_bases(l+1, p, tuple);
        }

        /* p -> p+1
           In this case we use p+1 symbols instead of p. Two cases occur. */
        tuple.p = p+1;
        /* Case 1: The new row in the parallel slit domain is inserted at the top, i.e. 
                   the transpositions 1, ..., l remain the same and we only insert the symbol 
                   p+1 in the (l+1)-th transposition, together with any symbol of 1, ..., p. */
        for(uint32_t i = p; i > 0; i--)
        {
            tuple[l+1] = Transposition(p+1, i);
            gen_bases(l+1, p+1, tuple);
        }

        /* Case 2: The new row in the parallel slit domain is not inserted at the top but at a position
                   i = 1, ..., p. Thus all indices i+1, ..., p are shifted up by one. */
        for( uint32_t i = p; i > 0; i-- )
        {
            Tuple tmp = tuple;
            for( uint32_t j = l; j > 0; j-- )
            {
                if( tmp[j].first >= i )
                {
                    tmp[j].first++;
                    if(tmp[j].second >= i )
                    {
                        tmp[j].second++;
                    }
                }
            }

            tmp[l+1] = Transposition(p+1, i);
            gen_bases(l+1, p+1, tmp);
        }

        /* p -> p+2
           Now we use p+2 symbols instead of p. Thus one row is inserted at the top of the 
           parallel slit domain, i.e. p+2 is the height of the (l+1)-th transposition. The other row is inserted 
           either directly below the top row or between the old rows. Since the symbol p+1 does not 
           occur in the transpositions 1, ..., l, both cases can be expressed by choosing a symbol
           i = 1, ..., p+1 and by shifting up all indices >= i by one. */ 
        tuple.p = p+2;
        for( uint32_t i = p+1; i > 0; i-- )
        {
            Tuple tmp = tuple;
            for( uint32_t j = l; j > 0; j-- )
            {
                if( tmp[j].first >= i )
                {
                    tmp[j].first++;
                    if(tmp[j].second >= i )
                    {
                        tmp[j].second++;
                    }
                }
            }

            tmp[l+1] = Transposition(p+2, i);
            gen_bases(l+1, p+2, tmp);
        }
    }
    else // Check whether the created h-tuple is really a generator, i.e. if it has the correct 
         // number of cycles. If this is the case, we add tuple to the basis elements of the 
         // p-th basis and store the index of tuple in this basis as the id of tuple.
    {
        if(tuple.permutation_type().num_cycles == m+1)
        {
            tuple.id = basis_complex[p].add_basis_element( tuple );
        }
    }
}

template< class MatrixComplex >
void MonoComplex< MatrixComplex > :: gen_differential(int32_t p)
{
    /**
     *  Instead of implementing the differential recursively, we use a direct formula to enumerate
     *  the sequences of indices in order to use threads.     
	 *  The sequences of indices we need to enumerate is given by the set 
	 *  \[ 
     *      \{(t_h, \ldots, t_1) \mid 0 \le t_q < q \,\}
	 *  \]
	 *  and for its enumeration we use that the map
	 *  \[ 
     *      \{0, \ldots, h! - 1\} = \{(t_h, \ldots, t_1) \mid 0 \le t_q < q \,\}
	 *  \]
	 *  \[
     *      k \mapsto \left( \left\lfloor \frac{k}{(q-1)!| \right\rfloor \pmod q\right)_q\,,
	 *  \]
	 *  is bijective. This is shown in the document s_qformel.pdf.	
    **/
    
    
    // Allocate enough space for the differential.
    // Todo: Test this.
    MatrixType differential( basis_complex[p-1].size(), basis_complex[p].size() );
    
    std::cout << "The map del o kappa_" << p << " is now computed..." << std::endl;
    // For each tuple t in the basis, we compute all basis elements that 
    // occur in kappa(t). 
    for( auto it : basis_complex[p].basis )
    {
        // parity of the exponent of the sign of the current summand of the differential
        int32_t parity = ((h*(h+1))/2) % 2;
        
        Tuple rand;
        uint32_t i;
        uint32_t s_q;
        
        for( uint32_t k = 0; k < factorial(h); k++ )
		// in each iteration we enumerate one sequence of indices according to the above formula        
		{
            Tuple neuer = it;
            bool norm_preserved = true;
            
            // Calculate phi_{(s_h, ..., s_1)}( Sigma )
            for( uint32_t q = 2; q <= h; q++ )
            {
                s_q = 1 + ( ( k / factorial(q-1)) % q );   
                parity += s_q % 2;
                if( neuer.phi(q, s_q) == false )
                {
                    norm_preserved = false;
                    break;
                }
            }
            
            // If phi_{(s_h, ..., s_1)}( Sigma ) is non-degenerate, we calculate the horizontal differential in .... and project back onto ....
            if( norm_preserved )   // Compute all horizontal boundaries.
            {
                for( i = 1; i < p; i++ )
                {
                    if( (rand = neuer.d_hor(i)) )
                    {
                        if( rand.monotone() == true ) // then it contributes to the differential with the computed parity
                        {
                            parity = (parity + i) % 2;
                            std::cout << it << "->" << rand << std::endl;
                            std::cout << it.id << "->" << basis_complex[p-1].id_of(rand) << " in " << "M_{" << basis_complex[p].size() << "," << basis_complex[p-1].size() << "} parity=" << parity << std::endl;
                            std::cout << std::endl;
//                            if (parity == 0)
//                            {
//                                differential(rand.id, it.id) += 1;
//                            }
//                            else
//                            {
//                                differential(rand.id, it.id) += -1;
//                            }
                        }
                    }
                }
            }
        }
    }
}

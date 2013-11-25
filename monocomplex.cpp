#include "monocomplex.hpp"

/*
 *
 *   MonoBasis
 *
 */

void MonoBasis :: add_basis_element (Tuple& t)
{
    basis.push_back(t);
}


uint64_t MonoBasis :: size()
{
    return basis.size();
}

/*
 *
 *   MonoComplex
 *
 */

MonoComplex :: MonoComplex(uint32_t _g, uint32_t _m) : g(_g), m(_m), h(2*_g + _m)
{
    Tuple tuple(h);

    tuple[1] = Transposition(2, 1);
    gen_bases(1, 2, tuple);  // We start with the transposition ... (2 1).
}

void MonoComplex :: gen_bases(uint32_t l, uint32_t p, Tuple& tuple)
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
    else // check whether the created h-tuple is really a generator, i.e. if it has the correct 
         // number of cycles
    {
        if(tuple.permutation_type().num_cycles == m+1)
        {
            basis_complex[p].add_basis_element( tuple );
        }
    }
}

void MonoComplex :: gen_differential(int32_t p)
{
    /**
        Instead of implementing the differential recursively, we use a direct formula to enumerate
        the sequences of indices in order to use threads. 
	The sequences of indices we need to enumerate is given by the set 
	\[ 
		\{(t_h, \ldots, t_1) \mid 0 \le t_q < q \,\}
	\]
	and for its enumeration we use that the map
	\[ 
		\{0, \ldots, h! - 1\} = \{(t_h, \ldots, t_1) \mid 0 \le t_q < q \,\}
	\]
	\[
		k \mapsto \left( \left\lfloor \frac{k}{(q-1)!| \right\rfloor \pmod q\right)_q\,,
	\]
	is bijective. This is shown in the document s_qformel.pdf.	
    **/

    uint64_t mono = 0;

    std::cout << "The map del o kappa_" << p << " is now computed...";
    #ifdef KAPPA_DEBUG_MONOKOMPLEX
         std::cout << std::endl;
    #endif

    // For each tuple t in the basis, we compute all basis elements that 
    // occur in kappa(t). 
    for( auto it : basis_complex[p] )
    {
        #ifdef KAPPA_DEBUG_MONOKOMPLEX
        std::cout << std::setfill('-') << std::setw(40) << "-" << std::setfill(' ') << std::endl
                  << "Consider tuple: " << it->to_string() << std::endl
                  << std::setfill('-') << std::setw(40) << "-" << std::setfill(' ') << std::endl;
        #endif

        // In order to use variables for each thread seperately with openmp, they have to be 
        // declared before the loop.
        Tuple neuer;
        Tuple rand;
        bool norm_preserved;
        uint32_t k;
        uint32_t i;
        uint32_t q;
        uint32_t s_q;
        #ifdef KAPPA_PARA
        #pragma omp parallel for private(neuer,rand,norm_preserved,k,i,q,s_q shared(mono)
        #endif
        for( k = 0; k < fakultaet(h); k++ )
	// in each iteration we enumerate one sequence of indices according to the above formula        
	{
            neuer = it;
            norm_preserved = true;

            #ifdef KAPPA_DEBUG_MONOKOMPLEX
            std::stringstream f_folge;
            std::string indexfolge = "1)";
            #endif

            // Calculate phi_{(s_h, ..., s_1)}( Sigma )
            for( q = 2; q <= h; q++ )
            {
                s_q = 1 + ( ( k / fakultaet(q-1)) % q );   

                #ifdef KAPPA_DEBUG_MONOKOMPLEX
                std::stringstream tmp;
                tmp << s_q;
                indexfolge = tmp.str() + "," + indexfolge;

                if( neuer.phi(q, s_q, &f_folge) == false )
                {
                    valid = false;
                    f_folge.str( std::string() );
                    break;
                }
                #else
                if( neuer.phi(q, s_q) == false )
                {
                    norm_preserved = false;
                    break;
                }
                #endif
            }
            
            // If phi_{(s_h, ..., s_1)}( Sigma ) is non-degenerate, we calculate the horizontal differential in .... and project back onto ....
            if( norm_preserved )   // Compute all horiyontal boundaries.
            {
                #ifdef KAPPA_DEBUG_MONOKOMPLEX
                indexfolge = "(" + indexfolge;
                std::cout << "Phi_" << indexfolge << "( " << it->to_string() << " ) = " << neuer.to_string() << std::endl;
                if( f_folge.str().length() > 1 )
                {
                    std::cout << "Thereby we have:" << std::endl
                              << f_folge.str();
                }
                f_folge.str( std::string() );
                bool bound_contr = false;
                #endif
                
                for( i = 1; i < p; i++ )
                {
                    if( (rand = neuer.d_hor(i)) )
                    {
                        #ifdef KAPPA_DEBUG_MONOKOMPLEX
                        bound_contr = true;
                        f_folge << "    The " << i << "-th boundary " << rand.to_string();
                        #endif
                        if( rand.monoton() == true ) // then it contributes to the differential
                        {
                            #ifdef KAPPA_PARA
			    // When we parallelize with openmp, mono++ must be an atomic operation.
                            #pragma omp atomic
                            #endif
                            mono++;
                            #ifdef KAPPA_DEBUG_MONOKOMPLEX
                            f_folge << " is monotone." << std::endl;
                            #endif
                        }
                        #ifdef KAPPA_DEBUG_MONOKOMPLEX
                        else
                        {
                            f_folge << " is not monotone." << std::endl;
                        }
                        #endif
                    }
                }
                #ifdef KAPPA_DEBUG_MONOKOMPLEX
                if( bound_contr == true )
                {
                    std::cout << "The horizontal boundaries are:" << std::endl
                              << f_folge.str() << std::endl;
                }
                else
                {
                    std::cout << "All horizontal boundaries are degenerate." << std::endl << std::endl;
                }
                #endif
            }
        }
    }

    #ifdef KAPPA_DEBUG_MONOKOMPLEX
    std::cout << "The number of monotonic pairs (t,t',i), where t' = del_i o kappa (t), is exactly " << mono << "." << std::endl;
    #else
    std::cout << mono << std::endl;
    #endif

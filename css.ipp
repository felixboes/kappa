#include "css.hpp"

template< class MatrixComplex >
ClusterSpectralSequence< MatrixComplex > :: ClusterSpectralSequence(uint32_t _g, uint32_t _m, SignConvention sgn) : g(_g), m(_m), h(2*_g + _m), sign_conv(sgn)
{
    Tuple tuple(h);
    tuple[1] = Transposition(2, 1);
    tuple.p = 2;
    
    gen_bases(1, 2, tuple);  // We start with the transposition ... (2 1).
}

template< class MatrixComplex >
void ClusterSpectralSequence< MatrixComplex > :: show_basis( int32_t p ) const
{
    if( basis_complex.count(p) )
    {
        std::cout << "This it the " << p << "-th basis:" << std::endl;
        std::cout <<  basis_complex.at(p) << std::endl;
    }
    else
    {
        std::cout << "The " << p << "-th basis is empty" << std::endl;
    }
}

template< class MatrixComplex >
void ClusterSpectralSequence< MatrixComplex > :: show_differential( int32_t p, uint32_t l ) const
{
    if( css_page.count(p) )
    {
        if( css_page[p].count(l) )
        {
            std::cout << "This it the " << l << "-th cluster of the " << p << "-th differential:" << std::endl;
            std::cout << css_page.at(p).at(l);
            std::cout << std::endl;
        }
    }
    else
    {
        std::cout << "The " << p << "-th differential is empty." << std::endl;
    }
}


template< class MatrixComplex >
void ClusterSpectralSequence< MatrixComplex > :: gen_bases(uint32_t s, uint32_t p, Tuple& tuple)
{
    /* Up to now we have determined all monotonic tuples of s transpositions containing the 
       symbols 1, ..., p, each at least once. We now add an (s+1)-th transposition and continue
       recursively.*/
    if(s < h) // There are h-s transpositionss+1eft to be determined.
    {
	/* From an s-tuple containing p symbols we can build up an (s+1)-tuple 
        with p, p+1 or p+2 symbols. */

        /* p -> p
           In this case we use the same number of symbols. Since we only
           enumerate monotonic tuples, the height of the (s+1)-th transposition needs to be p.
           We try out all possibilities for the second symbol in the (s+1)-th transposition. */
        tuple.p = p;
        for(uint32_t i = p-1; i > 0; i--)
        {
            tuple[s+1] = Transposition(p, i); 
            gen_bases(s+1, p, tuple);
        }

        /* p -> p+1
           In this case we use p+1 symbols instead of p. Two cases occur. */
        tuple.p = p+1;
        /* Case 1: The new row in the parallel slit domain is inserted at the top, i.e. 
                   the transpositions 1, ..., s remain the same and we only insert the symbol 
                   p+1 in the (s+1)-th transposition, together with any symbol of 1, ..., p. */
        for(uint32_t i = p; i > 0; i--)
        {
            tuple[s+1] = Transposition(p+1, i);
            gen_bases(s+1, p+1, tuple);
        }

        /* Case 2: The new row in the parallel slit domain is not inserted at the top but at a position
                   i = 1, ..., p. Thus all indices i+1, ..., p are shifted up by one. */
        for( uint32_t i = p; i > 0; i-- )
        {
            Tuple tmp = tuple;
            for( uint32_t j = s; j > 0; j-- )
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

            tmp[s+1] = Transposition(p+1, i);
            gen_bases(s+1, p+1, tmp);
        }

        /* p -> p+2
           Now we use p+2 symbols instead of p. Thus one row is inserted at the top of the 
           parallel slit domain, i.e. p+2 is the height of the (s+1)-th transposition. The other row is inserted 
           either directly below the top row or between the old rows. Since the symbol p+1 does not 
           occur in the transpositions 1, ..., s, both cases can be expressed by choosing a symbol
           i = 1, ..., p+1 and by shifting up all indices >= i by one. */ 
        tuple.p = p+2;
        for( uint32_t i = p+1; i > 0; i-- )
        {
            Tuple tmp = tuple;
            for( uint32_t j = s; j > 0; j-- )
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

            tmp[s+1] = Transposition(p+2, i);
            gen_bases(s+1, p+2, tmp);
        }
    }
    else // Check whether the created h-tuple is really a generator, i.e. if it has the correct 
         // number of cycles. If this is the case, we add tuple to the basis elements of the 
         // p-th basis and store the index of tuple in this basis as the id of tuple.
    {
        if(tuple.num_cycles() == m+1)
        {
            tuple.id = basis_complex[p].add_basis_element( tuple );
        }
    }
}

template< class MatrixComplex >
void ClusterSpectralSequence< MatrixComplex > :: gen_differentials( int32_t p )
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
    
    for( auto l : basis_complex.at(p) )
    {
        css_page[p][l] = MatrixType ( basis_complex[p-1].basis[l].size(), basis_complex[p].basis[l].size() );
        MatrixType& differential = css_page[p][l];
        // Initialize with zeros.
        differential.clear();
        
        // For each tuple t in the basis, we compute all basis elements that 
        // occur in kappa(t). 
        int32_t parity = 0;
        for( auto it : basis_complex[p].basis[l] )
        {
            Tuple boundary;
            uint32_t s_q;
            
            for( uint32_t k = 0; k < factorial(h); k++ )
            // in each iteration we enumerate one sequence of indices according to the above formula        
            {
                Tuple current_basis = it;
                bool norm_preserved = true;
                
                // parity of the exponent of the sign of the current summand of the differential
                if( sign_conv != no_signs ) 
                {
                    parity = ((h*(h+1))/2) % 2;
                }
                
                // Calculate phi_{(s_h, ..., s_1)}( Sigma )
                for( uint32_t q = 1; q <= h; q++ )
                {
                    s_q = 1 + ( ( k / factorial(q-1)) % q );   
                    if( sign_conv != no_signs )
                    {
                        parity += s_q;
                    }
                    if( current_basis.phi(q, s_q) == false )
                    {
                        norm_preserved = false;
                        break;
                    }
                }
                
                // If phi_{(s_h, ..., s_1)}( Sigma ) is non-degenerate, we calculate the horizontal differential in .... and project back onto ....
                if( norm_preserved )   // Compute all horizontal boundaries.
                {
                    std::map< uint8_t, int8_t > or_sign;
                    if( sign_conv == all_signs )
                    {
                        or_sign.operator =(std::move(current_basis.orientation_sign()));
                    }
    
                    for( uint32_t i = 1; i < p; i++ )
                    {
                        if( (boundary = current_basis.d_hor(i)) )
                        {
                            if( boundary.monotone() == true && boundary.num_cluster() == l ) // then it contributes to the differential with the computed parity
                            {
                                boundary.id = basis_complex[p-1].id_of(boundary);
                                
                                if( sign_conv == all_signs )
                                {
                                    int32_t actual_parity = (parity + i) % 2;
                                    if ( or_sign[i] == -1 )
                                    {
                                        actual_parity = (actual_parity + 1) % 2;
                                    }
                                    //std::cout << it << " " << i << ": The d^hor_i boundary of " << current_basis << ". This is " << boundary << std::endl;
                                    //std::cout << it.id << "->" << boundary.id << " in " << "M_{" << basis_complex[p-1].size() << "," << basis_complex[p].size() << "} parity=" << actual_parity << std::endl;
                                    //std::cout << std::endl;
                                    if ( actual_parity == 0 )
                                    {
                                        differential(boundary.id, it.id) += 1;
                                    }
                                    else
                                    {
                                        differential(boundary.id, it.id) += -1;
                                    }
                                }
                                else if( sign_conv == no_orientation_sign )
                                {
                                    if ( (parity + i) % 2 == 0 )
                                    {
                                        differential(boundary.id, it.id) += 1;
                                    }
                                    else
                                    {
                                        differential(boundary.id, it.id) += -1;
                                    }
                                }
                                else
                                {
                                    differential(boundary.id, it.id) += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template< class MatrixComplex >
void ClusterSpectralSequence< MatrixComplex >::erase_differential(int32_t p, uint32_t l)
{
    if( css_page.count(p) != 0 )
    {
        if( css_page[p].count(l) != 0)
        {
            // Delete matrix to save space.
            // Quote from the boost::ublas documentation http://www.boost.org/doc/libs/1_49_0/libs/numeric/ublas/doc/matrix.htm
            //
            // void resize (size_type size1, size_type size2, bool preserve = true)
            // Reallocates a matrix to hold size1 rows of size2 elements. The existing elements of the matrix are preseved when specified.
            css_page[p][l].resize(0,0,false);
            css_page[p].erase(l);
            if( css_page.count(p) == 0 )
            {
                css_page.erase(p);
            }
        }
    }
}

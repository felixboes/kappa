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
void MonoComplex< MatrixComplex > :: show_differential( int32_t p ) const
{
    if( matrix_complex.count(p) )
    {
        std::cout << "This it the " << p << "-th differential:" << std::endl;
        std::cout << matrix_complex.at(p);
        std::cout << std::endl;
    }
    else
    {
        std::cout << "The " << p << "-th differential is empty." << std::endl;
    }
}

template< class MatrixComplex >
void MonoComplex< MatrixComplex > :: show_differential_naive( int32_t p ) const
{
    if( matrix_complex_naive.count(p) )
    {
        std::cout << "This it the " << p << "-th differential (naive)." << std::endl;
        std::cout << matrix_complex_naive.at(p);
        std::cout << std::endl;
    }
    else
    {
        std::cout << "The " << p << "-th differential (naive) is empty." << std::endl;
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
    matrix_complex[p] = MatrixType ( basis_complex[p-1].size(), basis_complex[p].size() );
    MatrixType& differential = matrix_complex[p];
    
    std::cout << "The map pi o del o kappa_" << p << " is now computed..." << std::endl;
    // For each tuple t in the basis, we compute all basis elements that 
    // occur in kappa(t). 
    for( auto it : basis_complex[p].basis )
    {
        // parity of the exponent of the sign of the current summand of the differential
        int32_t parity = ((h*(h+1))/2) % 2;
        
        Tuple boundary;
        uint32_t s_q;
        
        for( uint32_t k = 0; k < factorial(h); k++ )
		// in each iteration we enumerate one sequence of indices according to the above formula        
		{
            Tuple current_basis = it;
            bool norm_preserved = true;
            parity = ((h*(h+1))/2) % 2;
            
            // Calculate phi_{(s_h, ..., s_1)}( Sigma )
            for( uint32_t q = 1; q <= h; q++ )
            {
                s_q = 1 + ( ( k / factorial(q-1)) % q );   
                parity += s_q;
                if( current_basis.phi(q, s_q) == false )
                {
                    norm_preserved = false;
                    break;
                }
            }
            
            // If phi_{(s_h, ..., s_1)}( Sigma ) is non-degenerate, we calculate the horizontal differential in .... and project back onto ....
            if( norm_preserved )   // Compute all horizontal boundaries.
            {
                for( uint32_t i = 1; i < p; i++ )
                {
                    if( (boundary = current_basis.d_hor(i)) )
                    {
                        if( boundary.monotone() == true ) // then it contributes to the differential with the computed parity
                        {
                            boundary.id = basis_complex[p-1].id_of(boundary);
                            
                            parity = (parity + i) % 2;
                            std::cout << it << " " << i << ": The d^hor_i boundary of " << current_basis << ". This is " << boundary << std::endl;
                            std::cout << it.id << "->" << boundary.id << " in " << "M_{" << basis_complex[p-1].size() << "," << basis_complex[p].size() << "} parity=" << parity << std::endl;
                            std::cout << std::endl;
                            if (parity == 0)
                            {
                                differential(boundary.id, it.id) += 1;
                            }
                            else
                            {
                                differential(boundary.id, it.id) += -1;
                            }
                        }
                    }
                }
            }
        }
    }
}


template< class MatrixComplex >
void MonoComplex< MatrixComplex > :: gen_differential_naive(int32_t p)
{    
    // Allocate enough space for the differential.
    // Todo: Test this.    
    matrix_complex_naive[p] = MatrixType ( basis_complex[p-1].size(), basis_complex[p].size() );
    
    std::cout << "The map pi o del o kappa_" << p << " (naive) is now computed..." << std::endl;
    // For each tuple t in the basis, we compute all basis elements that 
    // occur in kappa(t). 
    for( auto& it : basis_complex[p].basis )
    {
		//// here we recursively determine all contributing sequences s_h, ..., s_1.        

        std::vector<int32_t> s(h+1, 1);
		int pos = h;
		pi_del_phi_naive(it, s);
		
        while (pos > 0)
		{
			while( s[h] < h )
			{
				++s[h];
				pi_del_phi_naive(it, s);
			}
			while ( (pos > 0) && (s[pos] == pos) )
			{
				s[pos] = 1;
				--pos;
			}
			if (pos == 0)
            {
				break;
            }
            
			++s[pos];
			pi_del_phi_naive(it, s);
			pos = h;
		}
	}
}

// for one sequence s_p, ..., s_1, this calculates the multiple application of phi, of and of the d_i.
template< class MatrixComplex >
void MonoComplex< MatrixComplex > :: pi_del_phi_naive(const Tuple& it, std::vector<int32_t> & s)
{
    Tuple current_basis = it;
	Tuple boundary;
	bool norm_preserved = true;
    // parity of the exponent of the sign of the current summand of the differential
    int32_t parity = ((h*(h+1))/2) % 2;
    
    uint32_t& p = current_basis.p;
    MatrixType& differential = matrix_complex_naive[p];
    
	for (int32_t q = 1; q <= h; ++q)
	{
		parity += s[q] % 2;
		if (current_basis.phi(q, s[q]) == false)
		{
			norm_preserved = false;
			break;
		}
	}
    // If phi_{(s_h, ..., s_1)}( Sigma ) is non-degenerate, we calculate the horizontal differential in .... and project back onto ....
    if( norm_preserved )   // Compute all horizontal boundaries.
    {
        for( uint32_t i = 1; i < p; i++ )
        {
            if( (boundary = current_basis.d_hor_naive(i)) && boundary.monotone() ) // then it contributes to the differential with the computed parity
            {
                boundary.id = basis_complex[p-1].id_of(boundary);
                
                parity = (parity + i) % 2;
//                std::cout << it << "->" << boundary << std::endl;
//                std::cout << it.id << "->" << boundary.id << " in " << "M_{" << basis_complex[p-1].size() << "," << basis_complex[p].size() << "} parity=" << parity << std::endl;
//                std::cout << std::endl;
                if (parity == 0)
                {
                    differential(boundary.id, it.id) += 1;
                }
                else
                {
                    differential(boundary.id, it.id) += -1;
                }
            }
        }
    }
}


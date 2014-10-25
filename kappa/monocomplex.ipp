#include "monocomplex.hpp"

/*
 *
 *   MonoComplex
 *
 */
template< class MatrixComplex >
MonoComplex< MatrixComplex > :: MonoComplex(
        const uint32_t _g,
        const uint32_t _m,
        SignConvention sgn,
        const uint32_t number_working_threads,
        const uint32_t number_remaining_threads)
    : g(_g),
      m(_m),
      h(2*_g + _m),
      num_threads(number_working_threads + number_remaining_threads),
      sign_conv(sgn),
      matrix_complex(true)
{
    DiagonalizerType& diago = matrix_complex.get_diagonalizer();
    diago.num_working_threads = number_working_threads;
    diago.num_remaining_threads = number_remaining_threads;
    
    if ( Tuple::get_radial() ) // For radial cells, we have h = 2g + m - 1.
    {
        --h;
    }
    if (h == 0)
    {
        return;
    }
    // Generate all tuples with h transpositions containing the symbols 1, ..., p,
    // each at least once, with the correct number of cycles.
    Tuple tuple(h);
    tuple[1] = Transposition(2, 1);
    tuple.p = 2;
    gen_bases(1, 2, 1, tuple);  // We start with the transposition ... (2 1).
    // In the radial case, we also generate all tuples as above, but also containing
    // the symbol 0.
    if ( Tuple::get_radial() )
    {
        Tuple radial_tuple(h);
        radial_tuple.p = 1;
        radial_tuple[1] = Transposition(1, 0);
        gen_bases(1, 1, 0, radial_tuple);
    }
}

template< class MatrixComplex >
void MonoComplex< MatrixComplex > :: show_basis( const int32_t p ) const
{
    if( bases.count(p) )
    {
        std::cout << "This it the " << p << "-th basis:" << std::endl;
        const auto& basis_vector = bases.at(p).basis;
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
void MonoComplex< MatrixComplex > :: gen_bases( const uint32_t l, const uint32_t p, const uint32_t start_symbol, Tuple& tuple )
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
        for(uint32_t i = start_symbol; i < p; ++i)
        {
            tuple[l+1] = Transposition(p, i); 
            gen_bases(l+1, p, start_symbol, tuple);
        }

        /* p -> p+1
           In this case we use p+1 symbols instead of p. Two cases occur. */
        tuple.p = p+1;
        /* Case 1: The new row in the parallel slit domain is inserted at the top, i.e. 
                   the transpositions 1, ..., l remain the same and we only insert the symbol 
                   p+1 in the (l+1)-th transposition, together with any symbol of 1, ..., p. */
        for (uint32_t i = start_symbol; i <= p; ++i)
        {
            tuple[l+1] = Transposition(p+1, i);
            gen_bases(l+1, p+1, start_symbol, tuple);
        }

        /* Case 2: The new row in the parallel slit domain is not inserted at the top but at a position
                   i = 1, ..., p. Thus all indices i+1, ..., p are shifted up by one. */
        for(uint32_t i = start_symbol; i <= p; ++i)
        {
            Tuple tmp = tuple;
            for( uint32_t j = l; j >= 1; --j )
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

            tmp[l+1] = Transposition(p+1, i);
            gen_bases(l+1, p+1, start_symbol, tmp);
        }

        /* p -> p+2
           Now we use p+2 symbols instead of p. Thus one row is inserted at the top of the 
           parallel slit domain, i.e. p+2 is the height of the (l+1)-th transposition. The other row is inserted 
           either directly below the top row or between the old rows. Since the symbol p+1 does not 
           occur in the transpositions 1, ..., l, both cases can be expressed by choosing a symbol
           i = 1, ..., p+1 and by shifting up all indices >= i by one. */ 
        tuple.p = p+2;
        for( uint32_t i = start_symbol; i <= p + 1; ++i)
        {
            Tuple tmp = tuple;
            for( uint32_t j = l; j >= 1; j-- )
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

            tmp[l+1] = Transposition(p+2, i);
            gen_bases(l+1, p+2, start_symbol, tmp);
        }
    }
    else // Check whether the created h-tuple is really a generator, i.e. if it has the correct 
         // number of cycles. If this is the case, we add tuple to the basis elements of the 
         // p-th basis and store the index of tuple in this basis as the id of tuple.
    {
        if (tuple.has_correct_num_cycles(m))
        {
            tuple.id = bases[p].add_basis_element( tuple );
        }
    }
}


template <class MatrixType>
void update_differential(MatrixType &           differential,
                         const size_t           row,
                         const size_t           column,
                         const int32_t          parity,
                         const int8_t           i,
                         const int8_t           or_sign,
                         const SignConvention & sign_conv)
{
    differential(row, column) += typename MatrixType::CoefficientType( sign(parity, i, or_sign, sign_conv) );
}

template< class MatrixComplex >
void MonoComplex<MatrixComplex>::compute_boundary( Tuple & tuple, const uint32_t p, typename MatrixComplex::MatrixType & differential )
{
    int32_t parity = 0;
    Tuple boundary;
    uint32_t s_q;
    for( uint32_t k = 0; k < factorial(h); k++ )
    // in each iteration we enumerate one sequence of indices according to the above formula
    {
        Tuple current_basis = tuple;
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
                or_sign.operator =( std::move(current_basis.orientation_sign()) );
            }

            for( uint32_t i = Tuple::get_min_boundary_offset(); i <= p - Tuple::get_max_boundary_offset(); i++ )
            {
                if( (boundary = current_basis.d_hor(i)) )
                {
                    boundary.id = bases[p-1].id_of(boundary);
                    update_differential<MatrixType>(differential, tuple.id, boundary.id,
                                            parity, i, or_sign[i], sign_conv);
                }
            }
        }
    }
}

template< class MatrixComplex >
void monocomplex_work(
        MonoComplex<MatrixComplex> &            monocomplex,
        MonocomplexWork &                       work,
        const uint32_t                          p,
        typename MatrixComplex::MatrixType &    differential)
{
    for ( auto it : work)
    {
        monocomplex.compute_boundary( it, p, differential );
    }
}

template< class MatrixComplex >
void MonoComplex< MatrixComplex > :: gen_differential( const int32_t p )
{
    /**
     *  Instead of implementing the differential recursively, we use a direct formula to enumerate
     *  the sequences of indices in order to use threads.     
	 *  The sequences of indices we need to enumerate is given by the set 
	 *  \f[ 
     *      \{(t_h, \ldots, t_1) \mid 0 \le t_q < q \,\}
	 *  \f]
	 *  and for its enumeration we use that the map
	 *  \f[ 
     *      \{0, \ldots, h! - 1\} = \{(t_h, \ldots, t_1) \mid 0 \le t_q < q \,\}
	 *  \f]
	 *  \f[
     *      k \mapsto \left( \left\lfloor \frac{k}{(q-1)!}| \right\rfloor \pmod q\right)_q\,,
	 *  \f]
	 *  is bijective. This is shown in the document s_qformel.pdf.	
    **/

    // Allocate enough space for the differential.
    // Todo: Test this.
    MatrixType & differential = matrix_complex.get_current_differential();
    differential.resize( bases[p].size(), bases[p-1].size() );
    
    if( bases[p].size() == 0 || bases[p-1].size() == 0 )
    {
        return;
    }
    
    // For each tuple t in the basis, we compute all basis elements that
    // occur in kappa(t).
    std::vector<MonocomplexWork> elements_per_threads (num_threads);
    uint32_t num_elements_per_thread = bases[p].size() / num_threads;
    
    if (bases[p].size() % num_threads != 0)
    {
        ++num_elements_per_thread;
    }
    uint32_t t = 0;
    uint32_t cur = 0;
    for ( auto it : bases[p].basis )
    {
        elements_per_threads[t].push_back(it);
        ++cur;
        if ((t < num_threads - 1) and (cur == num_elements_per_thread * (t+1)))
        {
            ++t;
        }
    }
    std::vector<std::thread> workers(num_threads);
    for (uint32_t t = 0; t < num_threads; ++t)
    {
        workers[t] = std::thread(monocomplex_work<MatrixComplex>, std::ref(*this), std::ref(elements_per_threads[t]), p, std::ref(differential));
    }
    for (uint32_t t = 0; t < num_threads; ++t)
    {
        workers[t].join();
    }
}

template< class MatrixComplex >
typename MonoComplex< MatrixComplex >::MatrixType & MonoComplex< MatrixComplex > :: get_current_differential()
{
    return matrix_complex.get_current_differential();
}

template< class MatrixComplex >
const typename MonoComplex< MatrixComplex >::MatrixType & MonoComplex< MatrixComplex > :: get_current_differential() const
{
    return matrix_complex.get_current_differential();
}

template< class MatrixComplex >
void MonoComplex< MatrixComplex >::erase_current_differential()
{
    matrix_complex.erase();
}


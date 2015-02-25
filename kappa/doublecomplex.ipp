#include "doublecomplex.hpp"

template< class MatrixComplex >
DoubleComplex< MatrixComplex > :: DoubleComplex(
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
      matrix_complex(false)
{
    DiagonalizerType& diago = get_diagonalizer();
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
void DoubleComplex< MatrixComplex > :: show_basis( const int32_t p ) const
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
void DoubleComplex< MatrixComplex > :: gen_bases( const uint32_t l, const uint32_t p, const uint32_t start_symbol, Tuple& tuple )
{
    // Up to now we have determined all monotonic tuples of l transpositions containing the 
    // symbols 1, ..., p, each at least once. We now add an (l+1)-th transposition and continue
    // recursively.
    if(l < h) // There are h-l transpositions left to be determined.
    {
	// From an l-tuple containing p symbols we can build up an (l+1)-tuple 
    // with p, p+1 or p+2 symbols.

        // p -> p
        // In this case we use the same number of symbols.
        tuple.p = p;
        for(uint32_t a = start_symbol; a <= p; ++a)
        {
            for( uint32_t b = start_symbol; b < a; ++b )
            {
                tuple[l+1] = Transposition(a, b); 
                gen_bases(l+1, p, start_symbol, tuple);
            }
        }

        // p -> p+1
        // In this case we use p+1 symbols instead of p. Two cases occur.
        tuple.p = p+1;
        // Case 1: The new row in the parallel slit domain is inserted the spot a.
        for (uint32_t a = start_symbol; a <= p+1; ++a)
        {
            for( uint32_t b = start_symbol; b < a; ++b )
            {
                Tuple tmp = tuple;
                for( uint32_t j = l; j >= 1; j-- )
                {
                    if( tmp[j].first >= a )
                    {
                        tmp[j].first += 1;
                        if(tmp[j].second >= a )
                        {
                            tmp[j].second += 1;
                        }
                    }
                }
                tmp[l+1] = Transposition(a, b);
                gen_bases(l+1, p+1, start_symbol, tmp);
            }  
        }

        // Case 2: The new row in the parallel slit domain is inserted at the spot b.
        for (uint32_t a = start_symbol; a <= p+1; ++a)
        {
            for( uint32_t b = start_symbol; b < a; ++b )
            {
                Tuple tmp = tuple;
                for( uint32_t j = l; j >= 1; j-- )
                {
                    if( tmp[j].first >= b )
                    {
                        tmp[j].first += 1;
                        if(tmp[j].second >= b )
                        {
                            tmp[j].second += 1;
                        }
                    }
                }
                tmp[l+1] = Transposition(a, b);
                gen_bases(l+1, p+1, start_symbol, tmp);
            }  
        }

        // p -> p+2
        // Now we use p+2 symbols instead of p. Thus one row is inserted at the spot a and the other at the spot b.
        tuple.p = p+2;
        for (uint32_t a = start_symbol; a <= p+2; ++a)
        {
            for( uint32_t b = start_symbol; b < a; ++b )
            {
                Tuple tmp = tuple;
                for( uint32_t j = l; j >= 1; j-- )
                {
                    if( tmp[j].first >= b )
                    {
                        tmp[j].first += 1;
                        if( tmp[j].first >= a )
                        {
                            tmp[j].first +=1;
                        }
                        
                        if( tmp[j].second >= b )
                        {
                            tmp[j].second += 1;
                            if( tmp[j].second >= a )
                            {
                                tmp[j].second += 1;
                            }
                        }
                    }
                }
                tmp[l+1] = Transposition(a, b);
                gen_bases(l+1, p+2, start_symbol, tmp);
            }  
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
void update_differential_doublecomplex(MatrixType &           differential,
                         const size_t           row,
                         const size_t           column,
                         const int8_t           i,
                         const int8_t           or_sign,
                         const SignConvention & sign_conv)
{
    differential(row, column) += typename MatrixType::CoefficientType( sign(0, i, or_sign, sign_conv) );
}

template< class MatrixComplex >
void DoubleComplex<MatrixComplex>::compute_boundary( Tuple & tuple, const uint32_t p, typename MatrixComplex::MatrixType & differential )
{
    std::map< uint8_t, int8_t > or_sign;
    if( sign_conv == all_signs )
    {
        or_sign.operator =( std::move(tuple.orientation_sign()) );
    }

    for( uint32_t i = Tuple::get_min_boundary_offset(); i <= p - Tuple::get_max_boundary_offset(); i++ )
    {
        Tuple boundary = tuple.d_hor_double_complex(i);
        if( boundary )
        {
            boundary.id = bases[p-1].id_of(boundary);
            update_differential_doublecomplex<MatrixType>(differential, tuple.id, boundary.id, i, or_sign.at(i), sign_conv);
        }
    }
}

template< class MatrixComplex >
void doublecomplex_work(
        DoubleComplex<MatrixComplex> &          doublecomplex,
        DoubleComplexWork &                     work,
        const uint32_t                          p,
        typename MatrixComplex::MatrixType &    differential)
{
    for ( auto it : work)
    {
        doublecomplex.compute_boundary( it, p, differential );
    }
}

template< class MatrixComplex >
void DoubleComplex< MatrixComplex > :: gen_differential( const int32_t p )
{    
    MatrixType & differential = get_current_differential();  
    
    //differential.resize( bases[p].size(), bases[p-1].size() );
    differential.resize( ( bases.count(p) != 0 ? bases.at(p).size() : 0 ), ( bases.count(p-1) != 0 ? bases.at(p-1).size() : 0 ) );
    
    if( differential.size1() == 0 || differential.size2() == 0 )
    {
        return;
    }
    
    // For each tuple t in the basis, we compute all basis elements that
    // occur in kappa(t).
    std::vector<DoubleComplexWork> elements_per_threads (num_threads);
    uint32_t num_elements_per_thread = bases.at(p).size() / num_threads;
    
    if (bases.at(p).size() % num_threads != 0)
    {
        ++num_elements_per_thread;
    }
    uint32_t t = 0;
    uint32_t cur = 0;
    for ( auto it : bases.at(p).basis )
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
        workers[t] = std::thread(doublecomplex_work<MatrixComplex>, std::ref(*this), std::ref(elements_per_threads[t]), p, std::ref(differential));
    }
    for (uint32_t t = 0; t < num_threads; ++t)
    {
        workers[t].join();
    }
}

template< class MatrixComplex >
typename DoubleComplex< MatrixComplex >::MatrixType & DoubleComplex< MatrixComplex > :: get_current_differential()
{
    return matrix_complex.get_current_differential();
}

template< class MatrixComplex >
const typename DoubleComplex< MatrixComplex >::MatrixType & DoubleComplex< MatrixComplex > :: get_current_differential() const
{
    return matrix_complex.get_current_differential();
}

template< class MatrixComplex >
size_t DoubleComplex< MatrixComplex > :: num_rows() const
{
    return matrix_complex.num_rows();
}

template< class MatrixComplex >
size_t DoubleComplex< MatrixComplex > :: num_cols() const
{
    return matrix_complex.num_cols();
}

template< class MatrixComplex >
typename DoubleComplex< MatrixComplex >::DiagonalizerType & DoubleComplex< MatrixComplex > :: get_diagonalizer()
{
    return matrix_complex.get_diagonalizer();
}

template< class MatrixComplex >
const typename DoubleComplex< MatrixComplex >::DiagonalizerType & DoubleComplex< MatrixComplex > :: get_diagonalizer() const
{
    return matrix_complex.get_diagonalizer();
}

template< class MatrixComplex >
void DoubleComplex< MatrixComplex >::erase_current_differential()
{
    matrix_complex.erase();
}

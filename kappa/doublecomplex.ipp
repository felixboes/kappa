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
    
    if (h == 0)
    {
        return;
    }
    // Generate all highcells with h transpositions containing the symbols 1, ..., p,
    // each at least once, with the correct number of cycles.
    HighCell highcell(h);
    highcell[1] = Transposition(2, 1);
    highcell.p = 2;
    gen_bases(1, 2, 1, highcell);  // We start with the transposition ... (2 1).
}

template< class MatrixComplex >
void DoubleComplex< MatrixComplex > :: gen_bases( const uint32_t l, const uint32_t p, const uint32_t start_symbol, HighCell& highcell )
{
    // Up to now we have determined all monotonic highcells of l transpositions containing the 
    // symbols 1, ..., p, each at least once. We now add an (l+1)-th transposition and continue
    // recursively.
    if(l < h) // There are h-l transpositions left to be determined.
    {
	// From an l-highcell containing p symbols we can build up an (l+1)-highcell 
    // with p, p+1 or p+2 symbols.

        // p -> p
        // In this case we use the same number of symbols.
        highcell.p = p;
        for(uint32_t a = start_symbol; a <= p; ++a)
        {
            for( uint32_t b = start_symbol; b < a; ++b )
            {
                highcell[l+1] = Transposition(a, b); 
                gen_bases(l+1, p, start_symbol, highcell);
            }
        }

        // p -> p+1
        // In this case we use p+1 symbols instead of p. Two cases occur.
        highcell.p = p+1;
        // Case 1: The new row in the parallel slit domain is inserted the spot a.
        for (uint32_t a = start_symbol; a <= p+1; ++a)
        {
            for( uint32_t b = start_symbol; b < a; ++b )
            {
                HighCell tmp = highcell;
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
                HighCell tmp = highcell;
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
        highcell.p = p+2;
        for (uint32_t a = start_symbol; a <= p+2; ++a)
        {
            for( uint32_t b = start_symbol; b < a; ++b )
            {
                HighCell tmp = highcell;
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
    else // Check whether the created h-highcell is really a generator, i.e. if it has the correct 
         // number of cycles. If this is the case, we add highcell to the basis elements of the 
         // p-th basis and store the index of highcell in this basis as the id of highcell.
    {
        if (highcell.has_correct_num_cycles(m))
        {
            highcell.id = bases[p].add_basis_element( highcell );
            // find first occurence of monotony-break if any
            for( uint32_t i = 1; i <= h - 1; ++i )
            {
                if( highcell.at(i+1).first < highcell.at(i).first )
                {
                    HighCell boundary = highcell.d_ver(i);
                    boundary.id = bases[p].add_basis_element( boundary );
                }
            }
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
void DoubleComplex<MatrixComplex>::compute_boundary( HighCell & highcell, const uint32_t p, typename MatrixComplex::MatrixType & differential )
{
    std::map< uint8_t, int8_t > or_sign;
    if( sign_conv == all_signs )
    {
        or_sign.operator =( std::move(highcell.orientation_sign()) );
    }
    
    size_t offset = (bases.count(p) != 0 ? bases.at(p).size_red() : 0);
    for( uint32_t i = HighCell::get_min_boundary_offset(); i <= p - HighCell::get_max_boundary_offset(); i++ )
    {
        HighCell boundary = highcell.d_hor_double_complex(i);
        if( boundary && boundary.monotone() )
        {
            boundary.id = bases[p-1].id_of(boundary);
            update_differential_doublecomplex<MatrixType>(differential, highcell.id, offset + boundary.id, i, or_sign.at(i), sign_conv);
        }
    }
    
    for( int32_t j = 1; j < highcell.norm(); ++j )
    {
        HighCell boundary = highcell.d_ver(j);
        if( boundary )
        {
            boundary.id = bases[p].id_of(boundary);
            if( j % 2 == 0 )
            {
                differential(highcell.id, boundary.id) += -1;
            }
            else
            {
                differential(highcell.id, boundary.id) += 1;
            }
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
    differential.resize( ( bases.count(p) != 0 ? bases.at(p).size_col() : 0 ), ( bases.count(p) != 0 ? bases.at(p).size_red() : 0 ) + ( bases.count(p-1) != 0 ? bases.at(p-1).size_ess() : 0 ) );
    
    if( differential.size1() == 0 || differential.size2() == 0 )
    {
        return;
    }
    
    // For each highcell t in the basis, we compute all basis elements that
    // occur in kappa(t).
    std::vector<DoubleComplexWork> elements_per_threads (num_threads);
    uint32_t num_elements_per_thread = bases.at(p).size_col() / num_threads;
    
    if (bases.at(p).size_col() % num_threads != 0)
    {
        ++num_elements_per_thread;
    }
    uint32_t t = 0;
    uint32_t cur = 0;
    for ( auto it : bases.at(p).basis_col )
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

#include "css.hpp"

int32_t CSSBasis :: add_basis_element ( Tuple& t )
{
    uint32_t num_clusters = t.num_cluster();
    LBasisType& l_basis = basis[num_clusters];
    t.id = l_basis.size();
    l_basis.insert(t);
    
    return t.id;
}

int32_t CSSBasis :: size( const int32_t l ) const
{
    return (basis.count(l) != 0 ? basis.at(l).size() : 0 );
}

int32_t CSSBasis :: total_size() const
{
    int32_t size(0);
    for( auto& it : basis )
    {
        size += it.second.size();
    }
    return size;
}

int32_t CSSBasis :: id_of( Tuple& t ) const
{
    for( auto& l_basis_it : basis )
    {
        auto& l_basis = l_basis_it.second;
        auto it = l_basis.find(t);
        if( it != l_basis.end() )
        {
            return it->id;
        }
    }
    return -1;
}

int32_t CSSBasis :: total_id_of( Tuple& t ) const
{
    int64_t basis_offset(0);
    
    for( auto& l_basis_it : basis )
    {
        auto& l_basis = l_basis_it.second;
        auto it = l_basis.find(t);
        if( it != l_basis.end() )
        {
            return basis_offset + it->id;
        }
        basis_offset += this->size(l_basis_it.first);
    }
    return -1;
}

std::ostream& operator<< ( std::ostream& os, const CSSBasis& cb )
{
    for( auto& it : cb.basis )
    {
        os << "Cluster of size " << it.first << std::endl;
        for( auto& b : it.second )
        {
            os << b.id << ": " << b << std::endl;
        }
    }
    return os;
}


template< class MatrixComplex >
ClusterSpectralSequence< MatrixComplex > :: ClusterSpectralSequence(
        const uint32_t          genus,
        const uint32_t          num_punctures,
        const SignConvention    sgn,
        const uint32_t          number_working_threads,
        const uint32_t          number_remaining_threads)
    : g(genus), m(num_punctures), h(2*genus + num_punctures), num_threads(number_working_threads + number_remaining_threads), sign_conv(sgn), diff_complex(true)
{
    Tuple tuple(h);
    tuple[1] = Transposition(2, 1);
    tuple.p = 2;
    
    DiagonalizerType& diago = diff_complex.get_diagonalizer();
    diago.num_working_threads = number_working_threads;
    diago.num_remaining_threads = number_remaining_threads;
    
    gen_bases(1, 2, tuple);  // We start with the transposition ... (2 1).
}

template< class MatrixComplex >
void ClusterSpectralSequence< MatrixComplex > :: show_basis( const int32_t p ) const
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
void ClusterSpectralSequence< MatrixComplex > :: gen_bases( const uint32_t s, const uint32_t p, Tuple& tuple )
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
void ClusterSpectralSequence< MatrixComplex > :: gen_d0( const int32_t p, const int32_t l )
{
    MatrixType& differential = diff_complex.get_current_differential();
    differential.define_operations(MatrixType::main_and_secondary);
    differential.resize( basis_complex[p].basis[l].size(), basis_complex[p-1].basis[l].size(), true );
    differential.diagonal.clear();
    
    if( basis_complex[p].basis[l].size() == 0 || basis_complex[p-1].basis[l].size() == 0 )
    {
        return;
    }
    
    std::vector<CSSWork> elements_per_threads (num_threads);
    uint32_t num_elements_per_thread = basis_complex[p].basis[l].size() / num_threads;
    if (basis_complex[p].basis[l].size() % num_threads != 0)
    {
        ++num_elements_per_thread;
    }
    uint32_t t = 0;
    uint32_t cur = 0;
    for ( auto it : basis_complex[p].basis[l] )
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
        workers[t] = std::thread(css_work_0<MatrixComplex>, std::ref(*this), std::ref(elements_per_threads[t]), p, l, std::ref(differential));
    }
    for (uint32_t t = 0; t < num_threads; ++t)
    {
        workers[t].join();
    }
}

template< class MatrixComplex >
void css_work_0(ClusterSpectralSequence<MatrixComplex> & css,
                CSSWork & work,
                const int32_t p,
                const int32_t l,
                typename MatrixComplex::MatrixType & differential)
{
    for ( auto it : work)
    {
        css.gen_d0_boundary(it, p, l, differential);
    }
}

template< class MatrixComplex >
void ClusterSpectralSequence< MatrixComplex >::gen_d0_boundary(const Tuple & tuple,
                                                               const int32_t p,
                                                               const int32_t l,
                                                               typename MatrixComplex::MatrixType & differential)
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
                or_sign.operator=(std::move(current_basis.orientation_sign()));
            }

            for( int32_t i = 1; i < p; i++ )
            {
                if( (boundary = current_basis.d_hor(i)) )
                {
                    if( boundary.monotone() == true && boundary.num_cluster() == l ) // then it contributes to the differential with the computed parity
                    {
                        boundary.id = basis_complex[p-1].id_of(boundary);
                        update_differential(differential, tuple.id, boundary.id, parity, i, or_sign[i], sign_conv);
                    }
                }
            }
        }
    }
}

template< class MatrixComplex >
typename ClusterSpectralSequence< MatrixComplex >::MatrixType ClusterSpectralSequence< MatrixComplex >::gen_d1_row( const int32_t p, const int32_t l, const Tuple& basis_element)
{
    MatrixType single_row( 1, basis_complex[p-1].basis[l-1].size() ); // This will be a single row before applying row operations
    single_row.define_operations(MatrixType::OperationType::main_and_secondary);
    Tuple boundary;
    uint32_t s_q;
    
    if( basis_complex[p-1].basis[l-1].size() == 0 )
    {
        return single_row;
    }
    
    for( uint32_t k = 0; k < factorial(h); k++ )
    // in each iteration we enumerate one sequence of indices according to the above formula        
    {
        Tuple current_basis = basis_element;
        bool norm_preserved = true;
        
        int32_t parity = 0;
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
        
        // If phi_{(s_h, ..., s_1)}( Sigma ) is non-degenerate, we calculate the horizontal differential in
        //.... and project back onto ....
        if( norm_preserved )   // Compute all horizontal boundaries.
        {
            std::map< uint8_t, int8_t > or_sign;
            if( sign_conv == all_signs )
            {
                or_sign.operator=(std::move(current_basis.orientation_sign()));
            }

            for( int32_t i = 1; i < p; i++ )
            {
                if( (boundary = current_basis.d_hor(i)) )
                {
                    if( boundary.monotone() == true && boundary.num_cluster() == l - 1) // then it contributes to the differential with the computed parity
                    {
                        boundary.id = basis_complex[p-1].id_of(boundary);
                        update_differential(single_row, 0, boundary.id, parity, i, or_sign[i], sign_conv);
                    }
                }
            }
        }
    }
    return single_row;
}

template< class MatrixComplex >
void ClusterSpectralSequence< MatrixComplex >::gen_d1_apply_operations( MatrixType& row )
{
    MatrixType& differential = diff_complex.get_current_differential();
    const typename MatrixType::DiagonalType& diagonal = differential.diagonal;
    const size_t num_cols = row.size2();
    
    for( const auto& diag_pos : diagonal )
    {
        const auto& diag_row = diag_pos.first;
        const auto& diag_col = diag_pos.second;
        
        if( row.at( 0, diag_col ) != CoefficientType(0) )
        {
            CoefficientType lambda( - row(0, diag_col) / differential.main_at(diag_row, diag_col) );
            for( size_t j = diag_col; j < num_cols; ++j )
            {
                CoefficientType& a = row( 0, j );
                a += lambda * differential.main_at( diag_row, j );
            }
        }
    }
}

template< class MatrixComplex >
void css_work_1(ClusterSpectralSequence<MatrixComplex> & css,
              MonocomplexWork & work,
              const int32_t p,
              const int32_t l,
              typename MatrixComplex::MatrixType & differential,
              const size_t num_cols,
              const std::vector< size_t >& offset
              )
{
    for ( auto it : work)
    {
        typename MatrixComplex::MatrixType single_row = css.gen_d1_row( p, l, it );
        css.gen_d1_apply_operations( single_row );
        for( size_t j = 0; j < num_cols; ++j )
        {
            differential.sec_op( it.id, j ) = single_row( 0, j + offset[j] );
        }
    }
}

template< class MatrixComplex >
void ClusterSpectralSequence< MatrixComplex >::gen_d1_stage_1( const int32_t p, const int32_t l )
{
    MatrixType& differential = diff_complex.get_current_differential();
    const size_t num_rows = basis_complex[p].basis[l].size();
    const size_t num_cols = basis_complex[p-1].basis[l-1].size() - differential.diagonal.size();
    const size_t num_cols_diff = basis_complex[p-1].basis[l-1].size();
    differential.sec_resize(num_rows, num_cols, true );
    
    if( num_rows == 0 || num_cols == 0 )
    {
        return;
    }
    
    // The diagonal entries are sorted by column.
    typename MatrixType::DiagonalType diagonal_copy(differential.diagonal);
    
    // Compute offset
    std::vector< size_t > column_offset( num_cols, 0 );
    size_t offset = 0;
    
    size_t pos_diff = 0;
    auto diag_entry = diagonal_copy.cbegin();
    for( size_t pos_row = 0; pos_row < num_cols; ++pos_row, ++pos_diff )
    {
        // skip entries in the diagonal
        while( pos_diff < num_cols_diff && pos_diff == diag_entry->second )
        {
            ++pos_diff;
            ++diag_entry;
            ++offset;
        }
        
        // set pos_diff to the next position.
        column_offset[pos_row] = offset;
    }
    
    std::vector<CSSWork> elements_per_threads (num_threads);
    uint32_t num_elements_per_thread = basis_complex[p].basis[l].size() / num_threads;
    if (basis_complex[p].basis[l].size() % num_threads != 0)
    {
        ++num_elements_per_thread;
    }
    uint32_t t = 0;
    uint32_t cur = 0;
    for ( auto it : basis_complex[p].basis[l] )
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
        workers[t] = std::thread(css_work_1<MatrixComplex>, std::ref(*this), std::ref(elements_per_threads[t]), p, l, std::ref(differential), num_cols, std::cref(column_offset));
    }
    for (uint32_t t = 0; t < num_threads; ++t)
    {
        workers[t].join();
    }
}

template< class MatrixComplex >
void ClusterSpectralSequence< MatrixComplex >::prepare_d1_diag()
{
    diff_complex.get_current_differential().define_operations( MatrixType::OperationType::secondary );
}

template< class MatrixComplex >
void ClusterSpectralSequence< MatrixComplex >::erase_d0()
{
    MatrixType& differential = diff_complex.get_current_differential();
    differential.resize(0,0,true);
    differential.diagonal.clear();
}

template< class MatrixComplex >
void ClusterSpectralSequence< MatrixComplex >::erase_d1()
{
    MatrixType& differential = diff_complex.get_current_differential();
    differential.sec_resize(0,0,true);
}

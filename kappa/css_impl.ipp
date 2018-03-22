// The software kappa is a collection of programs to compute the homology of
// the moduli space of surfaces using the radial model.
// Copyright (C) 2013 - 2018  Felix Boes and Anna Hermann
// 
// This file is part of kappa.
// 
// kappa is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// kappa is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with kappa.  If not, see <http://www.gnu.org/licenses/>.


#include "css.hpp"

template< class MatrixComplex >
ClusterSpectralSequence< MatrixComplex > :: ClusterSpectralSequence(
        const uint32_t          genus,
        const uint32_t          num_punctures,
        const SignConvention    sgn,
        const uint32_t          number_working_threads,
        const uint32_t          number_remaining_threads)
    : g(genus),
      m(num_punctures),
      h(2*genus + num_punctures),
      num_threads(number_working_threads + number_remaining_threads),
      sign_conv(sgn),
      diff_complex(true)
{
    // Configure diagoanlizer
    DiagonalizerType& diago = diff_complex.get_diagonalizer();
    diago.num_working_threads = number_working_threads;
    diago.num_remaining_threads = number_remaining_threads;

    // Get parameter right.
    if (SymGrpTuple::is_radial() ) // For radial cells, we have h = 2g + m - 1.
    {
        --h;
    }
    if ( h == 0 )
    {
        return;
    }
    
    // Generate all tuples with h transpositions containing the symbols 1, ..., p,
    // each at least once, with the correct number of cycles.
    SymGrpTuple tuple(h);
    tuple[1] = Transposition(2, 1);
    tuple.p = 2;
    gen_bases(1, 2, 1, tuple);  // We start with the transposition ... (2 1).
    // In the radial case, we also generate all tuples as above, but also containing
    // the symbol 0.
    if (SymGrpTuple::is_radial() )
    {
        SymGrpTuple radial_tuple(h);
        radial_tuple.p = 1;
        radial_tuple[1] = Transposition(1, 0);
        gen_bases(1, 1, 0, radial_tuple);
    }
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
void ClusterSpectralSequence< MatrixComplex > :: gen_bases( const uint32_t s, const uint32_t p, const uint32_t start_symbol, SymGrpTuple& tuple )
{
    /* Up to now we have determined all monotonic tuples of s transpositions containing the 
       symbols 1, ..., p, each at least once. We now add an (s+1)-th transposition and continue
       recursively.*/
    if(s < h) // There are h-s transpositions left to be determined.
    {
	/* From an s-tuple containing p symbols we can build up an (s+1)-tuple 
        with p, p+1 or p+2 symbols. */

        /* p -> p
           In this case we use the same number of symbols. Since we only
           enumerate monotonic tuples, the height of the (s+1)-th transposition needs to be p.
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

template< class MatrixComplex >
void ClusterSpectralSequence< MatrixComplex > :: gen_d0( const int32_t p, const int32_t l )
{
    MatrixType& differential = diff_complex.get_current_differential();
    differential.define_operations(MatrixType::main_and_secondary);
    differential.resize( basis_complex[p].basis[l].size(), basis_complex[p-1].basis[l].size() );
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
void ClusterSpectralSequence< MatrixComplex >::gen_d0_boundary(const SymGrpTuple & tuple,
                                                               const int32_t p,
                                                               const int32_t l,
                                                               typename MatrixComplex::MatrixType & differential)
{
    int32_t parity = 0;
    SymGrpTuple boundary;
    uint32_t s_q;
    
    for( uint32_t k = 0; k < factorial(h); k++ )
    // in each iteration we enumerate one sequence of indices according to the above formula        
    {
        SymGrpTuple current_basis = tuple;
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

            for( uint32_t i = SymGrpTuple::get_min_boundary_offset(); i <= p - SymGrpTuple::get_max_boundary_offset(); i++ )
            {
                if( (boundary = current_basis.d_hor(i)) )
                {
                    if(boundary.fully_unstable() == true && boundary.num_clusters() == l ) // then it contributes to the differential with the computed parity
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
typename ClusterSpectralSequence< MatrixComplex >::MatrixType ClusterSpectralSequence< MatrixComplex >::gen_d1_row( const int32_t p, const int32_t l, const SymGrpTuple& basis_element)
{
    MatrixType single_row( 1, basis_complex[p-1].basis[l-1].size() ); // This will be a single row before applying row operations
    single_row.define_operations(MatrixType::OperationType::main_and_secondary);
    SymGrpTuple boundary;
    uint32_t s_q;
    
    if( basis_complex[p-1].basis[l-1].size() == 0 )
    {
        return single_row;
    }
    
    for( uint32_t k = 0; k < factorial(h); k++ )
    // in each iteration we enumerate one sequence of indices according to the above formula        
    {
        SymGrpTuple current_basis = basis_element;
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

            for( uint32_t i = SymGrpTuple::get_min_boundary_offset(); i <= p - SymGrpTuple::get_max_boundary_offset(); i++ )
            {
                if( (boundary = current_basis.d_hor(i)) )
                {
                    if(boundary.fully_unstable() == true && boundary.num_clusters() == l - 1) // then it contributes to the differential with the computed parity
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
              CSSWork & work,
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
            differential.sec_op( it.id, j ) = single_row.at( 0, j + offset[j] );
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
    differential.sec_resize(num_rows, num_cols );
    
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
        while( diag_entry != diagonal_copy.cend() && pos_diff < num_cols_diff && pos_diff == diag_entry->second )
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
    differential.resize(0,0);
    differential.diagonal.clear();
}

template< class MatrixComplex >
void ClusterSpectralSequence< MatrixComplex >::erase_d1()
{
    MatrixType& differential = diff_complex.get_current_differential();
    differential.sec_resize(0,0);
}

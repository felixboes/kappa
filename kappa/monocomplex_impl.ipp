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


#include "monocomplex.hpp"

/*
 *
 *   MonoComplex
 *
 */
template< class MatrixComplex >
MonoComplex< MatrixComplex > :: MonoComplex(
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
      diff_complex(false)
{
    // Configure diagoanlizer
    DiagonalizerType& diago = diff_complex.get_diagonalizer();
    diago.num_working_threads = number_working_threads;
    diago.num_remaining_threads = number_remaining_threads;

    // Get parameter right.
    if ( Tuple::get_radial() ) // For radial cells, we have h = 2g + m - 1.
    {
        --h;
    }
    if ( h == 0 )
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
void MonoComplex< MatrixComplex > :: gen_bases( const uint32_t s, const uint32_t p, const uint32_t start_symbol, Tuple& tuple )
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
            Tuple tmp = tuple;
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
            Tuple tmp = tuple;
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


template <class MatrixType>
void update_differential(MatrixType &           differential,
                         const size_t           row,
                         const size_t           column,
                         const int32_t          parity,
                         const int8_t           i,
                         const int8_t           or_sign,
                         const SignConvention & sign_conv)
{
    // differential(row, column) += typename MatrixType::CoefficientType( sign(parity, i, or_sign, sign_conv) );
    differential(row, column) += sign(parity, i, or_sign, sign_conv);
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
                    boundary.id = basis_complex[p-1].id_of(boundary);
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
    
    diff_complex.make_current_diff_old();
    MatrixType & differential = get_current_differential();  
    
    //differential.resize( basis_complex[p].size(), basis_complex[p-1].size() );
    differential.resize( ( basis_complex.count(p) != 0 ? basis_complex.at(p).size() : 0 ), ( basis_complex.count(p-1) != 0 ? basis_complex.at(p-1).size() : 0 ) );

    if( differential.size1() == 0 || differential.size2() == 0 )
    {
        return;
    }
    
    for( size_t j = 0; j < differential.size2(); ++j)
    {
        for( size_t i = 0; i < differential.size1(); ++i)
        {
            differential(i,j) = 0;
        }
    }

    // For each tuple t in the basis, we compute all basis elements that
    // occur in kappa(t).
    std::vector<MonocomplexWork> elements_per_threads (num_threads);
    uint32_t num_elements_per_thread = basis_complex.at(p).size() / num_threads;
    
    if (basis_complex.at(p).size() % num_threads != 0)
    {
        ++num_elements_per_thread;
    }
    uint32_t t = 0;
    uint32_t cur = 0;
    for ( auto it : basis_complex.at(p).basis )
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
void MonoComplex< MatrixComplex > :: apply_base_changes()
{
    diff_complex.apply_base_changes();
}

template< class MatrixComplex >
void MonoComplex< MatrixComplex > :: homchain(int32_t p, bool homology, int32_t maxdimension)
{
    if( homology == false )
    {
        if( basis_complex.empty() == true )
        {
            return;
        }

        if( homchain_cohomology_file.is_open() == false )
        {
            std::stringstream ss;
            ss << "homchain_cohomology_g_" << (int)g << "_m_" << (int)m << ".chomp";
            homchain_cohomology_file.open(ss.str());

            homchain_cohomology_file << "chain complex\n\nmaxdimension = " << (int)(1 + basis_complex.rbegin()->first) << "\n\n";
        }

        homchain_cohomology_file << "dimension " << p << "\n";
        auto mat = get_current_differential();
        size_t num_rows = mat.size1();
        size_t num_cols = mat.size2();
        CoefficientType coeff;

        for( size_t i = 0; i < num_rows; ++i)
        {
            coeff = 0;
            bool non_zero_column = false;
            homchain_cohomology_file << "   boundary a" << (int)(i + 1) << " = ";
            for( size_t j = 0; j < num_cols; ++j)
            {
                if( (bool)(coeff = mat(i,j)) )
                {
                    non_zero_column = true;
                    if( coeff > 0 )
                    {
                        homchain_cohomology_file << "+ " <<  coeff << " * a" << (int)(j + 1) << " ";
                    }
                    else
                    {
                        homchain_cohomology_file << "- " << -coeff << " * a" << (int)(j + 1) << " ";
                    }
                }
            }
            if( non_zero_column == false )
            {
                homchain_cohomology_file << "0";
            }
            homchain_cohomology_file << "\n";
        }
    }
    else
    {
        if( homchain_homology_file.is_open() == false )
        {
            std::stringstream ss;
            ss << "homchain_homology_g_" << (int)g << "_m_" << (int)m << ".chomp";
            homchain_homology_file.open(ss.str());

            homchain_homology_file << "chain complex\n\nmaxdimension = " << maxdimension << "\n\n";
        }

        homchain_homology_file << "dimension " << maxdimension - p << "\n";
        auto mat = get_current_differential();
        size_t num_rows = mat.size1();
        size_t num_cols = mat.size2();
        CoefficientType coeff;

        for( size_t j = 0; j < num_cols; ++j)
        {
            coeff = 0;
            bool non_zero_column = false;
            homchain_homology_file << "   boundary a" << (int)(j + 1) << " = ";
            for( size_t i = 0; i < num_rows; ++i)
            {
                if( (bool)(coeff = mat(i,j)) )
                {
                    non_zero_column = true;
                    if( coeff > 0 )
                    {
                        homchain_homology_file << "+ " <<  coeff << " * a" << (int)(i + 1) << " ";
                    }
                    else
                    {
                        homchain_homology_file << "- " << -coeff << " * a" << (int)(i + 1) << " ";
                    }
                }
            }
            if( non_zero_column == false )
            {
                homchain_homology_file << "0";
            }
            homchain_homology_file << "\n";
        }
    }
}


template< class MatrixComplex >
typename MonoComplex< MatrixComplex >::HomologyType MonoComplex< MatrixComplex > :: diagonalize_current_differential( const int32_t p, uint32_t max_possible_rank, const bool print_duration )
{
    if( max_possible_rank == 0 )
    {
        max_possible_rank = std::min( num_rows(), num_cols() );
    }
    atomic_uint& current_rank = get_diagonalizer().current_rank;
    
    // Compute the induced homology.
    atomic_uint state(0);   // Set state to 1 iff kernel and torsion are computed. This is done to terminate the 'monitoring thread'.
    ThisType* this_object = this;   // The C++ standard allows us to pass the this-pointer to the lambda function. Unfortunately, this may result a compiler bug....
    
    // Diagonalzing thread.
    auto partial_homology_thread = std::async( std::launch::async, [&]() -> HomologyType
    {
        // Always use one thread for diagonalizing at the moment!
        HomologyType ret = this_object->compute_current_kernel_and_torsion( p );
        state = 1;
        return ret;
    } );
    
    if( print_duration == true )
    {
        // Monitoring thread.
        auto monitor_thread = std::async( std::launch::async, [&]()
        {
            while( state != 1 )
            {
                std::cout << "Diagonalization " << current_rank << "/" << max_possible_rank << "\r";
                std::cout.flush();
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
        } );
        
        // Wait for threads to terminate.
        auto partial_homology = partial_homology_thread.get();
        monitor_thread.get();
        return partial_homology;
    }
    else
    {
        return partial_homology_thread.get();
    }
}

template< class MatrixComplex >
typename MonoComplex< MatrixComplex >::MatrixType & MonoComplex< MatrixComplex > :: get_current_differential()
{
    return diff_complex.get_current_differential();
}

template< class MatrixComplex >
const typename MonoComplex< MatrixComplex >::MatrixType & MonoComplex< MatrixComplex > :: get_current_differential() const
{
    return diff_complex.get_current_differential();
}

template< class MatrixComplex >
size_t MonoComplex< MatrixComplex > :: num_rows() const
{
    return diff_complex.num_rows();
}

template< class MatrixComplex >
size_t MonoComplex< MatrixComplex > :: num_cols() const
{
    return diff_complex.num_cols();
}

template< class MatrixComplex >
typename MonoComplex< MatrixComplex >::DiagonalizerType & MonoComplex< MatrixComplex > :: get_diagonalizer()
{
    return diff_complex.get_diagonalizer();
}

template< class MatrixComplex >
const typename MonoComplex< MatrixComplex >::DiagonalizerType & MonoComplex< MatrixComplex > :: get_diagonalizer() const
{
    return diff_complex.get_diagonalizer();
}

template< class MatrixComplex >
void MonoComplex< MatrixComplex >::erase_current_differential()
{
    diff_complex.erase();
}

template< class MatrixComplex >
typename MonoComplex< MatrixComplex >::HomologyType MonoComplex< MatrixComplex > :: compute_current_kernel_and_torsion( const int32_t n )
{
    return diff_complex.compute_current_kernel_and_torsion(n);
}

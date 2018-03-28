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


#include "sym_grp_tuple.hpp"
#include "sym_grp_ehr_bases_generator.hpp"

bool SymGrpTuple::radial = false;
uint32_t SymGrpTuple::min_symbol = 1;
uint32_t SymGrpTuple::min_boundary_offset = 1;
uint32_t SymGrpTuple::max_boundary_offset = 1;

SymGrpTuple :: SymGrpTuple() :
    p(0),
    rep()
{
}

SymGrpTuple :: SymGrpTuple(const size_t num_entries) :
    p(0),
    rep( num_entries, Transposition(0, 0) )
{
}

SymGrpTuple :: SymGrpTuple(const uint32_t symbols, const size_t num_entries) :
    p(symbols),
    rep( num_entries, Transposition(0, 0) )
{
}

void SymGrpTuple::parallel_case()
{
    radial = false;
    min_symbol = 1;
    min_boundary_offset = 1;
    max_boundary_offset = 1;
}

void SymGrpTuple::radial_case()
{
    radial = true;
    min_symbol = 0;
    min_boundary_offset = 0;
    max_boundary_offset = 0;
}

bool SymGrpTuple::is_radial()
{
    return radial;
}

uint32_t SymGrpTuple::get_min_symbol()
{
    return min_symbol;
}

uint32_t SymGrpTuple::get_min_boundary_offset()
{
    return min_boundary_offset;
}

uint32_t SymGrpTuple::get_max_boundary_offset()
{
    return max_boundary_offset;
}

Transposition & SymGrpTuple :: at(const size_t n)
{
    return rep.at(n-1);
}

const Transposition & SymGrpTuple :: at(const size_t n) const
{
    return rep.at(n-1);
}

Transposition& SymGrpTuple :: operator[](const size_t n)
{
    return rep[n-1];
}

int32_t SymGrpTuple :: num_entries() const
{
    return rep.size();
}

bool SymGrpTuple :: operator==(const SymGrpTuple& t) const
{
    // Observe how two vectors are compared
    // Operations == and != are performed by first comparing sizes, and if they match, the elements are compared sequentially
    // using algorithm equal, which stops at the first mismatch.
    // Source: http://www.cplusplus.com/reference/vector/vector/operators/
    return (this->rep == t.rep);
}

bool SymGrpTuple :: operator!=(const SymGrpTuple& t) const
{
    // Observe how two vectors are compared
    // Operations == and != are performed by first comparing sizes, and if they match, the elements are compared sequentially
    // using algorithm equal, which stops at the first mismatch.
    // Source: http://www.cplusplus.com/reference/vector/vector/operators/
    return (this->rep != t.rep);
}

SymGrpTuple :: operator bool() const
{
    if (rep.size() == 0)
    {
        return false;
    }
    return true;
}

bool SymGrpTuple :: has_correct_num_cycles(const size_t m) const
{
    // min_symbol = 1 iff radial = false iff num_components = m + 'exactly one boundary component of the surface'
    return num_cycles() == m + min_symbol;
}

bool SymGrpTuple :: is_multiple_of_a() const
{
    if( this->operator bool() == true )
    {
        if( this->at( this->rep.size() ) == Transposition( this->p, this->p-1 ) )
        {
            const size_t rep_size_min_1 = rep.size() - 1;
            const size_t p_min_2 = this -> p - 2;
            for( size_t i = 0; i < rep_size_min_1; ++i )
            {
                if( rep.at(i).first > p_min_2 || rep.at(i).second > p_min_2 )
                {
                    return false;
                }
            }
            return true;
        }
    }

    return false;
}

uint32_t SymGrpTuple :: num_cycles() const
{
    // count instead the number of cycles of sigma_out_inv
    return PermutationManager::num_cycles(sigma_out_inv());
}

SymGrpTuple::ConnectedComponents SymGrpTuple::connected_components() const
{
    // TODO Adapt for radial case
    // Compare http://www.boost.org/doc/libs/1_49_0/libs/graph/example/connected_components.cpp
    typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;

    // Build graph.
    Graph G;
    for( auto edge : rep )
    {
        boost::add_edge(edge.first, edge.second, G);
    }

    // Let boost compute the connected compontents.
    ConnectedComponents components(p+1);
    boost::connected_components(G, &components[0]);

    return components;
}

/*
std::vector< size_t > SymGrpTuple::shuffle_positions() const
{
    std::vector< size_t > slits;
    const size_t h = num_entries();
    size_t k = 2*h;
    slits.push_back(k);
    
    while( k != 0 )
    {
        const size_t j = (k+1)/2;
        const size_t a = ( 2*j - k == 1 ? at(j).first : at(j).second );
        
        std::cout << j << " -> " << a << std::endl;
        
        // search for slits on the same hight but below a (i.e. longer slits).
        for( size_t l = j-1; l >= 1; --l )
        {
            if( a == at(l).first )
            {
                k = 2*l;
                // [n3290: 6.1/1]: [..] The scope of a label is the function in which it appears. [..]
                // Therefore we can use the name 'end_outer_loop'.
                goto end_outer_loop;
            }
            else if( a == at(l).second )
            {
                k = 2*l-1;
                // [n3290: 6.1/1]: [..] The scope of a label is the function in which it appears. [..]
                // Therefore we can use the name 'end_outer_loop'.
                goto end_outer_loop;
            }
        }
        
        // if a is the longest slit at the first line, we are done.
        if( a == 1 )
        {
            k = 0;
            // [n3290: 6.1/1]: [..] The scope of a label is the function in which it appears. [..]
            // Therefore we can use the name 'end_outer_loop'.
            goto end_outer_loop;
        }
        
        // find the first slit below a.
        for( size_t l = h; l >= 1; --l )
        {
            if( a-1 == at(l).first )
            {
                k = 2*l;
                // [n3290: 6.1/1]: [..] The scope of a label is the function in which it appears. [..]
                // Therefore we can use the name 'end_outer_loop'.
                goto end_outer_loop;
            }
            else if( a-1 == at(l).second )
            {
                k = 2*l-1;
                // [n3290: 6.1/1]: [..] The scope of a label is the function in which it appears. [..]
                // Therefore we can use the name 'end_outer_loop'.
                goto end_outer_loop;
            }
        }
        
        end_outer_loop:
        slits.push_back(k);
    }
    return slits;
}

std::vector< size_t > SymGrpTuple::slits() const
{
    std::vector< size_t > the_slits;
    size_t h = num_entries();
    
    for( size_t i = 1; i <= p; ++i )
    {
        for( size_t j = 1; j <= h; ++j )
        {
            if( at(j).first == i )
            {
                the_slits.push_back( 2*j );
            }
            else if ( at(j).second == i )
            {
                the_slits.push_back( 2*j - 1 );
            }
        }
    }
    
    return the_slits;
}

SymGrpTuple SymGrpTuple::Q_term( const std::vector< size_t >& shuffle_slit_conf, const size_t num_shuffle_pos, const size_t fuse_pos ) const
{
    const size_t h = num_entries();
    SymGrpTuple res(2*p, 2*num_entries()); // We reduce p later.
    const auto the_slits = slits();
    const auto the_shuffles = shuffle_positions();
    
    for( size_t j = 0; j < 2*h; ++j )
    {
                                                                         
    }
    
    return res;
}
*/

int32_t SymGrpTuple::num_clusters() const
{
    typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;

    // Build graph.
    Graph G;
    for( auto edge : rep )
    {
        boost::add_edge(edge.first, edge.second, G);
    }

    // Let boost compute the connected compontents.
    ConnectedComponents components(p+1);
    // returns the number of components.
    // Here 0 is a valid vertex, so we must shrink the number of components if radial = false.
    return boost::connected_components(G, &components[0]) - min_symbol;
}

bool SymGrpTuple :: fully_unstable() const
{
    // A SymGrpTuple is fully unstable iff the sequence of all at(i) is monotone.
    for( int32_t i = 1; i <= num_entries() - 1; i++ )
    {
        if( at(i+1).first < at(i).first )
        {
            return false;
        }
    }
    return true;
}

bool SymGrpTuple :: f(const uint32_t i)
{
    // Denote g_{i+1} | g_i by (ab)(cd).
    uint8_t a = at(i+1).first;
    uint8_t b = at(i+1).second;
    uint8_t c = at(i).first;
    uint8_t d = at(i).second;

    if( a == c ) // (a )(a )
    {
        if( b == d ) // (ab)(ab) = id and (a*)(a*) = 0
        {
            rep = std::vector<Transposition>();
            return false;
        }
        else if( a == b ) // (a*)(ad) = (ad*) -> (d*)(ad)
        {
            at(i+1) = Transposition(d,d);
            return true;
        }
        else if( a == d ) // (ab)(a*) = (ab*) -> (b*)(ab)
        {
            at(i+1) = Transposition(b,b);
            at(i)   = Transposition(a,b);
            return true;
        }
        else // (ab)(ad) = (adb) = (bd)(ab) or (db)(ab)
        {
            at(i+1) = Transposition(std::max(b,d), std::min(b,d));
            at(i)   = Transposition(a,b);
            return true;
        }
    }
    else if( a == b ) // (a*)(  )
    {
        if( a == d ) // (a*)(ca) = (ca*) -> (a*)(ca)
        {
            return true;
        }
        else // (a*)(c*) or (a*)(cd) and the two transpositions commute
        {
            if( a < c )
            {
                return true;
            }
            else
            {
                std::swap( at(i+1), at(i));
                return true;
            }
        }
    }
    // Case: a > b
    else if( b == c ) // (ab)(b )
    {
        if( c == d ) // (ab)(b*) = (ab*) -> (b*)(ab)
        {
            std::swap( at(i+1), at(i) );
            return true;
        }
        else // (ab)(bd) = (abd) and (abd)' = (ad) and (abd)(ad) = (bd)
        {
            at(i+1) = Transposition(b,d);
            at(i)   = Transposition(a,d);
            return true;
        }
        // These are all cases since a > b >= d.
    }
    else if( a == d ) // (ab)(ca) = (acb) and c > a -> (ab)(ca)
    {
        return true;
    }
    else if( b == d ) // (ab)(cb) = (abc) -> a < c: (ab)(cb) oder a > c:(cb)(ac)
    {
        if( a < c )
        {
            return true;
        }
        else    // (abc)' = (ac) and (abc)(ac) = (cb)
        {
            at(i+1) = Transposition(c,b);
            at(i)   = Transposition(a,c);
            return true;
        }
    }
    else // (ab)(c ) with d not a and not b, hence disjoint transpositions
    {
        if( a < c )
        {
            return true;
        }
        else
        {
            std::swap( at(i+1), at(i) );
            return true;
        }
    }
    std::cerr << "Error in 'bool SymGrpTuple::f(uint32_t i)' -> Reached impossible case." << std::endl;
    return false;
}
                 
bool SymGrpTuple :: phi( const uint32_t q, const uint32_t i )
{
    if( i == 0 || i > q )
    {
        std::cerr << "Error in 'bool SymGrpTuple::phi( uint32_t q, uint32_t i)' -> q=" << q << " and i=" << i << std::endl;
        return false;
    }

    for( uint32_t j = q-1; j >= i; j-- ) // The loop terminates due to i > 0.
    {
        if( f(j) == false ) // The norm of the product falls.
        {
            return false;
        }
    }
    return true;
}

SymGrpTuple SymGrpTuple :: d_hor( const uint8_t k ) const
{
    SymGrpTuple boundary = *this;
    
    // start with sigma_0.
    Permutation sigma = long_cycle();
    Permutation sigma_inv = long_cycle_inv();
    
    for(uint8_t q = 1; q <= boundary.num_entries(); ++q)
    {
        // Write tau_q = (a,b)
        auto a = boundary[q].first;
        auto b = boundary[q].second;
        
        //Write (k, sigma_{q-1}(k)) = (k,l)
        auto l = sigma[k];
        
         // Compute tau':
        // Most of the time the transpositions are disjoint hence (a,b)(k,l) = (k,l)(a,b) and
        // the left transposition will be killed by D_k
        if( k != a && k != b && l != a && l != b )
        {
        }
        // The degenerate case:
        // k and l are both part of tau_q.
        else if( a == std::max(k, l) and ( b == std::min(k, l) or k == l))
        {
            return SymGrpTuple();
        }
        // The non degenerate case:
        else
        {   
            // Compute Z = (a,b)(k,l)
            // Z(k) = l iff l != a and l != b hence k == a or k == b
            // In this case: (Z(k),k)Z = (k,l)(a,b)(k,l) = (c,l) with c != k
            if( l != a && l != b)
            {
                if(a != k)
                {
                    boundary[q].first  = std::max(a,l);
                    boundary[q].second = std::min(a,l);
                }
                else
                {
                    boundary[q].first  = std::max(b,l);
                    boundary[q].second = std::min(b,l);
                }
            }
            // Z(k) != l iff l = a or l = b hence k != a and k != b
            // In this case: (Z(k),k)Z = (c,k)(a,b)(k,l) = (c,l) with c != l,
            // but we see that this is just (a,b):
            // if l != a, then k,l != a and we map
            // a to a to b != c since otherwise a would map to k
            // therefore a = c and l = b.
            // Thus we do not need to alter the boundary.
        }
        
        // Compute sigma_{q}
        // (a,b)sigma(k) =
        //   a          if k = sigma^{-1}(b)
        //   b          if k = sigma^{-1}(a)
        //   sigma(k)   else
        // this is done by swapping the values of sigma^{-1}(a) and sigma^{-1}(b) under sigma:
        std::swap( sigma[ sigma_inv[a] ], sigma[ sigma_inv[b] ] );
        
        // sigma^{-1}(a,b) (k) = 
        //   sigma^{-1}(b)  if k = a
        //   sigma^{-1}(a)  if k = b
        //   sigma^{-1}(k)  else
        // this is done by interchanging the values of a and b under sigma^{-1}
        std::swap( sigma_inv[a], sigma_inv[b] );
    }

    // boundary is fully_unstable iff its renormalization is fully_unstable. Thus we check for unstability now to avoid
    // unnecessary renormalization.
    if (not boundary.fully_unstable())
    {
        return SymGrpTuple();
    }
    
    // Renormalize all tau'
    for(uint8_t q = 1; q <= boundary.num_entries(); ++q)
    {
        // Write tau_q = (a,b)
        auto& a = boundary[q].first;
        auto& b = boundary[q].second;

        if(a > k)
        {
            a--;
        }
        if(b > k)
        {
            b--;
        }
    }
    
    boundary.p -= 1;
    
    return boundary;
}

SymGrpTuple SymGrpTuple :: d_hor_reduced(const uint8_t i) const
{
    const SymGrpTuple boundary = this->d_hor(i);
    if( boundary.is_multiple_of_a() == false )
    {
        return boundary;
    }
    else
    {
        return SymGrpTuple();
    }
}

SymGrpTuple SymGrpTuple :: d_hor_double_complex( const uint8_t k ) const
{
    SymGrpTuple boundary = *this;
    
    // start with sigma_0.
    Permutation sigma = long_cycle();
    Permutation sigma_inv = long_cycle_inv();
    
    for(uint8_t q = 1; q <= boundary.num_entries(); ++q)
    {
        // Write tau_q = (a,b)
        auto a = boundary[q].first;
        auto b = boundary[q].second;
        
        //Write (k, sigma_{q-1}(k)) = (k,l)
        auto l = sigma[k];
        
         // Compute tau':
        // Most of the time the transpositions are disjoint hence (a,b)(k,l) = (k,l)(a,b) and
        // the left transposition will be killed by D_k
        if( k != a && k != b && l != a && l != b )
        {
        }
        // The degenerate case:
        // k and l are both part of tau_q.
        else if( a == std::max(k, l) and ( b == std::min(k, l) or k == l))
        {
            return SymGrpTuple();
        }
        // The non degenerate case:
        else
        {   
            // Compute Z = (a,b)(k,l)
            // Z(k) = l iff l != a and l != b hence k == a or k == b
            // In this case: (Z(k),k)Z = (k,l)(a,b)(k,l) = (c,l) with c != k
            if( l != a && l != b)
            {
                if(a != k)
                {
                    boundary[q].first  = std::max(a,l);
                    boundary[q].second = std::min(a,l);
                }
                else
                {
                    boundary[q].first  = std::max(b,l);
                    boundary[q].second = std::min(b,l);
                }
            }
            // Z(k) != l iff l = a or l = b hence k != a and k != b
            // In this case: (Z(k),k)Z = (c,k)(a,b)(k,l) = (c,l) with c != l,
            // but we see that this is just (a,b):
            // if l != a, then k,l != a and we map
            // a to a to b != c since otherwise a would map to k
            // therefore a = c and l = b.
            // Thus we do not need to alter the boundary.
        }
        
        // Compute sigma_{q}
        // (a,b)sigma(k) =
        //   a          if k = sigma^{-1}(b)
        //   b          if k = sigma^{-1}(a)
        //   sigma(k)   else
        // this is done by swapping the values of sigma^{-1}(a) and sigma^{-1}(b) under sigma:
        std::swap( sigma[ sigma_inv[a] ], sigma[ sigma_inv[b] ] );
        
        // sigma^{-1}(a,b) (k) = 
        //   sigma^{-1}(b)  if k = a
        //   sigma^{-1}(a)  if k = b
        //   sigma^{-1}(k)  else
        // this is done by interchanging the values of a and b under sigma^{-1}
        std::swap( sigma_inv[a], sigma_inv[b] );
    }
    
    // Renormalize all tau'
    for(uint8_t q = 1; q <= boundary.num_entries(); ++q)
    {
        // Write tau_q = (a,b)
        auto& a = boundary[q].first;
        auto& b = boundary[q].second;

        if(a > k)
        {
            a--;
        }
        if(b > k)
        {
            b--;
        }
    }
    
    boundary.p -= 1;
    
    return boundary;
}

Permutation SymGrpTuple::sigma_out_inv() const
{
    Permutation sigma_inv = long_cycle_inv();

    for (uint8_t i = 1; i <= num_entries(); ++i) {
        PermutationManager::multiply(sigma_inv, at(i));
    }

    return sigma_inv;
}

Permutation SymGrpTuple::sigma_out() const
{
    return PermutationManager::inverse(sigma_out_inv());
}

std::map< uint8_t, int8_t > SymGrpTuple::orientation_sign( ) const
{
    Permutation sigma = sigma_out();
    return PermutationManager::orientation_sign_for_ehrenfried(sigma);
}

Permutation SymGrpTuple::long_cycle() const
{
    return Permutation::long_cycle(p);
}

Permutation SymGrpTuple::long_cycle_inv() const
{
    return Permutation::long_cycle_inv(p);
}

SymGrpTuple create_cell( const size_t num_entries, ... )
{
    SymGrpTuple t(num_entries);
    t.p = 0;
    va_list args;
    va_start(args, num_entries);
    
    for ( size_t i = num_entries; i > 0; --i )
    {
        const uint8_t a = va_arg(args, int);
        const uint8_t b = va_arg(args, int);
        t[i] = Transposition( a, b );
        t.p = std::max( t.p, (uint32_t)std::max( a, b ) );
    }
    va_end(args);
    return t;
}

size_t HashSymGrpTuple :: operator ()( const SymGrpTuple &tuple ) const
{
    size_t hashvalue = 0;
    size_t offset = 2;
    for( const auto& cit : tuple.rep )
    {
        hashvalue += offset*(cit.first + 8*cit.second);
        offset *= 16;
    }
    return hashvalue;
}

std::ostream& operator<< (std::ostream& stream, const SymGrpTuple& tuple)
{
    for( auto it = tuple.rep.crbegin(); it != tuple.rep.crend(); ++it )
    {
        stream << "(" << (uint32_t)it->first << ",";
        if( it->first != it->second )
        {
            stream << (uint32_t)it->second << ")";
        }
        else
        {
            stream << "*)";
        }
    }
    
    return stream;
}

SymGrpTuple operator*( const SymGrpTuple& v_2, const SymGrpTuple& v_1 )
{
    const uint32_t & p_1 = v_1.p;
    const uint32_t & p_2 = v_2.p;
    
    const size_t h_1 = v_1.num_entries();
    const size_t h_2 = v_2.num_entries();
    const size_t h_prod = h_1 + h_2;
    SymGrpTuple prod;
    prod.p = p_1 + p_2;
    prod.id = std::numeric_limits< size_t >::max();
    
    auto& rep = prod.rep;
    rep.reserve( h_prod * sizeof( Transposition ) );
    
    for( size_t i = 0; i < h_1; ++i )
    {
        rep.emplace_back( v_1.rep.at(i) );
    }
    for( size_t i = 0; i < h_2; ++i )
    {
        rep.emplace_back( Transposition( p_1 + v_2.rep.at(i).first, p_1 + v_2.rep.at(i).second ) );
    }
    
    return prod;
}

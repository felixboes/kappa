#include "tuple.hpp"

Permutation::Permutation()
    : data()
{
    // intentionally do nothing
}

Permutation::Permutation(uint8_t size, uint8_t init)
    : data(size, init)
{
    // intentionally do nothing
}

Permutation::Permutation(const Permutation & other)
    :
      data(other())
{
    // intentionally do nothing
}

std::vector<uint8_t> Permutation::operator()() const
{
    return data;
}

uint8_t & Permutation::operator[](uint8_t i)
{
    return data[i];
}

uint8_t const & Permutation::at(uint8_t i) const
{
    return data.at(i);
}

uint8_t Permutation::size() const
{
    return data.size();
}

bool Permutation::is_contained(uint8_t i) const
{
    return (data[i] != 0);
}

bool Permutation::is_fix_point(uint8_t i) const
{
    return (data[i] == i);
}

std::ostream& operator<< (std::ostream& stream, const Permutation& permutation)
{
    stream << "permutation: " << std::endl;
    for (uint8_t i = 0; i < permutation.size(); ++i)
    {
        if (permutation.is_contained(i))
        {
            stream << (size_t) i << " maps to " << permutation.at(i) << std::endl;
        }
    }
    return stream;
}

Tuple :: Tuple() :
    p(0),
    rep()
{
}

Tuple :: Tuple(size_t h) :
    p(0),
    rep( h, Transposition(0, 0) )
{
}

Tuple :: Tuple(uint32_t symbols, size_t h) :
    p(symbols),
    rep( h, Transposition(0, 0) )
{
}

Transposition& Tuple :: at(size_t n)
{
    return rep[n-1];
}

Transposition const & Tuple :: at(size_t n) const
{
    return rep[n-1];
}

Transposition& Tuple :: operator[](size_t n)
{
    return at(n);
}

int32_t Tuple :: norm() const
{
    if( this->operator bool() == true )
    {
        return rep.size();
    }
    else // The Tuple is not valid.
    {
        return 0;
    }
}

bool Tuple :: operator==(const Tuple& t) const
{
    // Observe how two vectors are compared
    // Operations == and != are performed by first comparing sizes, and if they match, the elements are compared sequentially
    // using algorithm equal, which stops at the first mismatch.
    // Source: http://www.cplusplus.com/reference/vector/vector/operators/
    return (this->rep == t.rep);
}

bool Tuple :: operator!=(const Tuple& t) const
{
    // Observe how two vectors are compared
    // Operations == and != are performed by first comparing sizes, and if they match, the elements are compared sequentially
    // using algorithm equal, which stops at the first mismatch.
    // Source: http://www.cplusplus.com/reference/vector/vector/operators/
    return (this->rep != t.rep);
}

Tuple :: operator bool() const
{
    if (rep.size() == 0)
    {
        return false;
    }
    return true;
}

std::ostream& operator<< (std::ostream& stream, const Tuple& tuple)
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
    
    std::cout << " Basis number (if any): " << tuple.id;
    
    return stream;
}

uint32_t Tuple :: num_cycles(size_t min_symbol)
{
    // Since the t_i are transpositions, one can instead count the number of cycles of (p p-1 ... 1) t_1 ... t_h.
    uint32_t num_cycles = 0;
    Permutation sigma_inv = long_cycle_inv();

     // multiply with t_1, ..., t_h
    for( int32_t i = 1; i <= norm(); i++ )
    {
        uint8_t k = at(i).first;
        uint8_t l = at(i).second;
        std::swap( sigma_inv[ k ], sigma_inv[ l ] );
    }

    // count the cycles
    std::vector<bool> visited( p+1, false );    // visited[0] is not used
    for( uint8_t i = min_symbol; i <= p; ) // We iterate through all cycles and mark the used symbols.
    {
        // consider the next cycle
        num_cycles += 1;
        visited[i] = true;
        uint8_t j = sigma_inv[i];
        
        while( j != i ) // mark all symbols in this cycle
        {
            visited[j] = true;
            j = sigma_inv[j];
        }

        // find the next unvisited cycle
        for( ++i; i <= p && visited[i]; ++i )
        {
        }
    }
    return num_cycles;
}

Tuple::ConnectedComponents Tuple::connected_components() const
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
    int32_t num = boost::connected_components(G, &components[0]);
    components[0] = num - 1;

    return components;
}

int32_t Tuple::num_cluster() const
{
    return connected_components()[0];
}

bool Tuple :: monotone()
{
    // A tuple is monotone iff the sequence of all at(i) is monotone.
    for( int32_t i = 1; i <= norm() - 1; i++ )
    {
        if( at(i+1).first < at(i).first )
        {
            return false;
        }
    }
    return true;
}

bool Tuple :: f(uint32_t i)
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
            std::swap( at(i), at(i) );
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
    std::cerr << "Error in 'bool Tuple::f(uint32_t i)' -> Reached impossible case." << std::endl;
    return false;
}
                 
bool Tuple :: phi( uint32_t q, uint32_t i )
{
    if( i == 0 || i > q )
    {
        std::cerr << "Error in 'bool Tuple::phi( uint32_t q, uint32_t i)' -> q=" << q << " and i=" << i << std::endl;
        return false;
    }

    for( uint32_t j = q-1; j >= i; j-- ) // The loop terminates due to i > 0.
    {
        #ifdef KAPPA_DEBUG_TUPEL
        std::cerr << "    f_" << j << "( " << *this << " ) = ";
        #endif
        if( f(j) == false ) // The norm of the product falls.
        {
            #ifdef KAPPA_DEBUG_TUPEL
            std::cerr << *this << std::endl;
            #endif
            return false;
        }
        #ifdef KAPPA_DEBUG_TUPEL
        else
        {
            std::cerr << *this << std::endl;
        }
        #endif
    }
    return true;
}

Tuple Tuple :: d_hor( uint8_t k ) const
{
    Tuple boundary = *this;
    
    // start with sigma_0.
    Permutation sigma = long_cycle();
    Permutation sigma_inv = long_cycle_inv();
    
    for(uint8_t q = 1; q <= boundary.norm(); ++q)
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
            return Tuple();
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

    // boundary is monotone iff its renormalization is monotone. Thus we check for monotony now to avoid unnecessary renormalization.
    if (not boundary.monotone())
    {
        return Tuple();
    }
    
    // Renormalize all tau'
    for(uint8_t q = 1; q <= boundary.norm(); ++q)
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
    return boundary;
}

Tuple Tuple :: d_hor_naive( uint8_t i ) const
{
    if( i == 0 || i >= p )
    {
        return Tuple();
    }
    Tuple boundary = *this;
    Permutation sigma = long_cycle();            // sigma = sigma_0
    Permutation sigma_inv = long_cycle_inv();    // sigma^{-1}
    
    // Compute the permutations sigma_i. The only symbols that change when we multiply with tau_i = (k,l) are
    // sigma_{i-1}^{-1}(k) and sigma_{i-1}^{-1}(l).
    for( int32_t l = 1; l <= boundary.norm(); l++ )
    {
        // write tau_l = (a,b)
        uint8_t a = boundary.at(l).first;
        uint8_t b = boundary.at(l).second;
        
        // Compute D_i(sigma_{l-1}^{-1})
        Permutation sigma_inv_tmp(sigma_inv);
        sigma_inv_tmp[ sigma[i] ] = sigma_inv[i];
        sigma_inv_tmp[i] = i;
        
        // compute sigma_l = tau_l sigma_{l-1}
        sigma[ sigma_inv[a] ] = b;
        sigma[ sigma_inv[b] ] = a;
        
        // compute sigma_l^{-1}
        std::swap( sigma_inv[ a ], sigma_inv[ b ] );
        
        // Compute D_i(sigma_l)
        Permutation sigma_tmp(sigma);
        sigma_tmp[ sigma_inv[i] ] = sigma[i];
        sigma_tmp[i] = i;
        
        // D_i(sigma_l) is degenerate iff sigma_l(i) = i
        if( sigma[i] == i )
        {
            return Tuple();
        }
        
        for( uint8_t j = 1; j <=p; ++j )
        {
            // This case occures at least once, otherwise tau_l = id.
            if( sigma_tmp[ sigma_inv_tmp[j] ] != j )
            {
                // Compute the new transpositon tau''_l = (c,d).
                uint8_t c = j;
                uint8_t d = sigma_tmp[ sigma_inv_tmp[j] ];
                
                if( c > i )
                {
                    c--;
                }
                if( d > i )
                {
                    d--;
                }
                if(c == d)
                {
                    std::cerr << "Error in 'Tuple Tuple :: d_hor_naiv(uint8_t i)' -> Reached impossible case." << std::endl;
                }
                boundary.at(l) = Transposition ( std::max(c,d),  std::min(c,d) );
            }
        }
    }

    return boundary;
}

Permutation Tuple::sigma_h() const
{
    // initialize with sigma_0
    Permutation sigma_inv = long_cycle_inv();

    for (uint8_t i = 1; i <= norm(); ++i)
    {
        // write tau_i = (a, b)
        const uint8_t& a = at(i).first;
        const uint8_t& b = at(i).second;
        
        // compute sigma_i^{-1}
        std::swap( sigma_inv[ a ], sigma_inv[ b ] );
    }
    
    // compute sigma
    Permutation sigma(p+1, 0);
    for( uint8_t i = 0; i <= p; ++i )
    {
        sigma[sigma_inv[i]] = i;
    }
    
    return sigma;
}

/**
 * Determines the decompositions of sigma into cycles including fix points.
 * \return A map consisting of all cycles of pi with their smallest element as key (except for
 * the liniel, for which the key is p)
 */
std::map< uint8_t, Permutation > Tuple::cycle_decomposition ( const Permutation & pi ) const
{
    std::map<uint8_t, Permutation> cycle_decomp;
    std::vector<bool> visited(p+1, false);
    for( uint8_t i = 1; i <= p; ) // We iterate through all cycles and mark the used symbols.
    {
        // determine the cycle of i.
        Permutation cycle(p+1, 0);
        
        uint8_t prev;     // previous symbol
        uint8_t cur = i; // current symbol
        bool liniel_found = false;
        
        do // mark all symbols in this cycle
        {
            visited[cur] = true;
            prev = cur;
            cur = pi.at(prev);
            cycle[prev] = cur;
            if ( prev == p )
            {
                liniel_found = true;
            }
        }while( cur != i );
        // note that since the for-loop runs ascendingly, the smallest element of the cycle
        // is i
        if ( liniel_found == true )
        {
            cycle_decomp[p] = cycle;
        }
        else
        {
            cycle_decomp[i] = cycle;
        }
        // find the next unvisited cycle
        for( ++i; i <= p && visited[i]; ++i )
        {
        }
    }
    return cycle_decomp;
}

std::map< uint8_t, int8_t > Tuple::orientation_sign( ) const
{
    Permutation sigma = sigma_h();

    std::map< uint8_t, Permutation > cycle_decomp = cycle_decomposition(sigma);

    std::map< uint8_t, int8_t > sign;
    // set the sign to 1 for all elements of the cycle of p
    for ( size_t i = 1; i < cycle_decomp.at(p).size(); ++i)
    {
        if (cycle_decomp.at(p)[i] != 0) // then i belongs to this cycle
        {
            sign[i] = 1;
        }
    }
    
    uint8_t i = 1; // counter of cycles
    for ( auto it_1 = cycle_decomp.begin(); it_1 != cycle_decomp.end(); ++it_1 )
    {
        // The liniel is treated seperately and is stored in cycle_decomp[p]
        if( it_1->first == p )
        {
            continue;
        }
    
        Permutation cycle = it_1->second;
        uint8_t min_symbol = it_1->first;
        // if the cycle is a fixpoint (a), we set sign(a) = 0 for the sake of completeness.
        if (cycle.is_fix_point(min_symbol))
        {
            sign[min_symbol] = 0;
            ++i;
            continue;
        }
        uint8_t second_min_symbol = 0;
        // determine the second min symbol of the cycle
        for (size_t m = 1; m < cycle.size(); ++m)
        {
            if (cycle[m] != 0 && m != min_symbol)
            {
                second_min_symbol = m;
                break;
            }
        }
        // Find k.
        // note that
        //   a_{1,1} < ... < a_{i,1} < b ,
        // hence
        //   a_{i-l,1} < b    for    l >= 0
        // is impossible and we can start to search at the position k = i.
        uint8_t k = i;
        
        // note that since we exclude the case that second_min_sybols will be sorted in at the end,
        // min_symbol can be sorted in between the cycles and this loop will never reach the liniel.
        
        // note that initially (for k = i) it_2.first = a_{k,1} < b and we found our position k iff the first time b < a_{k+1,1}
        // in the other case we have again a_{k+1,1} < b. Induction.
        // Moreover the case 'b > a_{m,1}' is also included as the liniel is stored at cycle_decomp[p] and b < p.
        for ( auto it_2 = it_1; it_2 != cycle_decomp.end(); ++it_2 )
        {
            uint8_t next_min_symbol = (std::next(it_2))->first;
            if ( second_min_symbol < next_min_symbol )
            {
                if ( ((k - i) % 2) == 0 )
                {
                    sign[min_symbol] = 1;
                }
                else
                {
                    sign[min_symbol] = -1;
                }
                break;
            }
            ++k;
        }
        
        // for all other symbols of the cycle, we set the sign to 1
        for (size_t c = 1; c < cycle.size(); ++c)
        {
            if (cycle[c] != 0 and c != min_symbol)
            {
                sign[c] = 1;
            }
        }
        ++i;
    }
    return sign;
}

Permutation Tuple::long_cycle() const
{
    Permutation sigma(p+1, 0);
    for(uint8_t k = 0; k < p; ++k)
    {
        sigma[k] = k+1;
    }
    sigma[p] = 0;
    return sigma;
}

Permutation Tuple::long_cycle_inv() const
{
    Permutation sigma(p+1, 0);
    for(uint8_t k = 1; k <= p; ++k)
    {
        sigma[k] = k-1;
    }
    sigma[0] = p;
    return sigma;
}

size_t HashTuple :: operator ()( const Tuple &tuple ) const
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


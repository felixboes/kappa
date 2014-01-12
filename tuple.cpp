#include "tuple.hpp"



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
    return this->rep == t.rep;
}

bool Tuple :: operator!=(const Tuple& t) const
{
    // Observe how two vectors are compared
    // Operations == and != are performed by first comparing sizes, and if they match, the elements are compared sequentially
    // using algorithm equal, which stops at the first mismatch.
    // Source: http://www.cplusplus.com/reference/vector/vector/operators/
    return this->rep != t.rep;
}

Tuple :: operator bool() const
{
    for( auto const &it : rep )
    {
        if( it.first == 0 || it.second == 0 )
        {
            return false;
        }
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

PermutationType Tuple :: permutation_type()
{
    // Since the t_i are transpositions, one can instead count the number of cycles of (p p-1 ... 1) t_1 ... t_h.
    PermutationType pt;   
    Tuple::Permutation sigma_inv = long_cycle_inv();

    // multiply with t_1, ..., t_h and count punctures
    for( uint32_t i = 1; i <= norm(); i++ )
    {
        if(at(i).first == at(i).second)
        {
            pt.num_punctures += 1;
        }

        uint8_t k = at(i).first;
        uint8_t l = at(i).second;
        std::swap( sigma_inv[ k ], sigma_inv[ l ] );
    }

    // count the cycles
    std::vector<bool> visited( p+1, false );    // visited[0] is not used
    for( uint8_t i = 1; i <= p; ) // We iterate through all cycles and mark the used symbols.
    {
        // consider the next cycle
        pt.num_cycles += 1;
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
    return pt;
}

bool Tuple :: monotone()
{
    // A tuple is monotone iff the sequence of all at(i) is monotone.
    for( uint32_t i = 1; i <= norm() - 1; i++ )
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
    // Compare Mehner: page 49, computed by a clear case distinction.
    if( i == 0 || i >= norm() ) // i has to be smaller than h.
    {
        std::cerr << "Error in 'bool Tuple :: f(uint32_t i)' -> i=" << i << " >= h=" << norm() << std::endl;
        return false;
    }
    else    // Denote g_{i+1} | g_i by (ab)(cd).
    {
        uint8_t a = at(i+1).first;
        uint8_t b = at(i+1).second;
        uint8_t c = at(i).first;
        uint8_t d = at(i).second;
        
        if( a == c ) // (a )(a )
        {
            if( b == d ) // (ab)(ab) = id and (a*)(a*) = 0
            {
                rep.assign( norm(), Transposition(0,0) );
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
    if(0 == k || k >= p)
    {
        return Tuple(this->norm());
    }
    
    Tuple boundary = *this;
    
    // start with sigma_0.
    Tuple::Permutation sigma = long_cycle();
    Tuple::Permutation sigma_inv = long_cycle_inv();
    
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
        // The degenerate cases:
        // k and k = sigma_{q-1}(k) are both part of tau_q.
        else if( k == l && ( a == k || b == k ) )
        {
            return Tuple();
        }
        // tau_q = (k, sigma_{q-1}(k)).
        else if( a == std::max(k, l) && b == std::min(k, l) )
        {
            return Tuple(boundary.norm());
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
    for(uint8_t q = 1; q <= boundary.norm(); ++q)
    {
        // Write tau_q = (a,b)
        auto& a = boundary[q].first;
        auto& b = boundary[q].second;
        
        if(a == k || b == k)
        {
            std::cerr << "Error in 'Tuple Tuple::d_hor( uint32_t k )' -> Reached impossible case." << std::endl;
        }
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
        return Tuple(this->norm());
    }
    Tuple boundary = *this;
    Tuple::Permutation sigma = long_cycle();            // sigma = sigma_0
    Tuple::Permutation sigma_inv = long_cycle_inv();    // sigma^{-1}
    
    // Compute the permutations sigma_i. The only symbols that change when we multiply with tau_i = (k,l) are
    // sigma_{i-1}^{-1}(k) and sigma_{i-1}^{-1}(l).
    for( uint32_t l = 1; l <= boundary.norm(); l++ )
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
            return Tuple( boundary.norm() );
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

/**
 * Determines the decompositions of sigma into cycles including fix points.
 * \return A map consisting of all cycles of sigma with their smallest element as key (except for
 * the liniel, for which the key is p)
 * \warning Since Mehner's example in the orientation chapter looks as if fix points were considered
 * in the cycle decomposition, we also consider them here.
 * \todo We assume that the liniel is not a fix point. Is that always true?
 */
std::map< uint8_t, Tuple::Permutation > Tuple::cycle_decomposition ( const Tuple::Permutation& sigma ) const
{
    std::map<uint8_t, Tuple::Permutation> cycles;
    std::vector<bool> visited(p+1, false);
    for( uint8_t i = 1; i <= p; ) // We iterate through all cycles and mark the used symbols.
    {
        // consider the next cycle
        Tuple::Permutation cycle;
        visited[i] = true;
        uint8_t j = i;
        uint8_t k = sigma.at(i);
        cycle[j] = k;
        bool liniel_found = false;
        while( k != i ) // mark all symbols in this cycle
        {
            visited[k] = true;
            j = k;
            k = sigma.at(k);
            cycle[j] = k;
            if ( j == p)
            {
                liniel_found = true;
            }
        }
        // note that since the for-loop runs ascendingly, the smallest element of the cycle
        // is i
        if ( liniel_found == true )
        {
            cycles[p] = cycle;
        }
        else
        {
            cycles[i] = cycle;
        }
        // find the next unvisited cycle
        for( ++i; i <= p && visited[i]; ++i )
        {
        }
    }
    return cycles;
}

std::map< uint8_t, int8_t > Tuple::orientation_sign( const Tuple::Permutation& sigma ) const
{
    std::map< uint8_t, Tuple::Permutation > cycles = cycle_decomposition(sigma);
    std::map< uint8_t, int8_t > sign;
    // set the sign to 1 for all elements of the cycle of p
    for ( auto &it : cycles.at(p) )
    {
        uint8_t k = it.first;
        sign[k] = 1;
    }

    uint8_t i = 1; // counter of cycles
    for ( auto &it_1 = cycles.begin(); it_1.first != p; ++it_1 )
    {
        Permutation cycle = it_1.second;
        
        uint8_t min_symbol = it_1.first;
        
        // if the cycle is a fixpoint (a), we set sign(a) = 0.
        if ( cycle.size() == 1)
        {
            sign[min_symbol] = 0;
        }
        else 
        {
            // for the minimum symbol of the cycle, we set the sign according to the formula
            uint8_t second_min_symbol = (std::next(cycle.begin()))->first;

            // \Note The liniel does not seem to appear in Mehner's cycle decomposition, thus
            // we need to ingore it
            if ( second_min_symbol > (std::prev(std::prev(cycle.end()))->first) )
            {
                uint8_t k = cycles.size() - 1;
                if ( k - i % 2)
                {
                    sign[min_symbol] = -1;
                }
                else
                {
                    sign[min_symbol] = 1;
                }
            }
            else
            {
                uint8_t k = 1;
                // note that since we exclude the case that second_min_sybols will be sorted in at the end,
                // min_symbol can be sorted in between the cycles and this loop will never reach the liniel.
                for ( auto &it_2 : cycles)
                {
                    uint8_t left = it_2.first;
                    uint8_t right = (std::next(&it_2))->first;
                    if ( left < min_symbol && min_symbol < right)
                    {
                        if ( k - i % 2)
                        {
                            sign[min_symbol] = -1;
                        }
                        else
                        {
                            sign[min_symbol] = 1;
                        }
                        break;
                    }
                    ++k;
                }
            }
            // for all other symbols of the cycle, we set the sign to 1
            for ( auto &it_2 = ++cycle.begin(); it_2 != cycle.end(); ++it_2)
            {
                uint8_t k = it_2->first;
                sign[k] = 1;
            }
        }
        ++i;
    }
}

Tuple::Permutation Tuple::long_cycle() const
{
    Tuple::Permutation sigma;
    for(uint8_t k = 1; k < p; ++k)
    {
        sigma[k] = k+1;
    }
    sigma[p] = 1;
    return sigma;
}

Tuple::Permutation Tuple::long_cycle_inv() const
{
    Tuple::Permutation sigma;
    for(uint8_t k = 2; k <= p; ++k)
    {
        sigma[k] = k-1;
    }
    sigma[1] = p;
    return sigma;
}

void Tuple :: print_permutation(Permutation sigma) const
{
    std::vector<bool> visited(p+1, false);
    for( uint8_t i = 1; i <= p; ) // We iterate through all cycles and mark the used symbols.
    {
        // consider the next cycle
        visited[i] = true;
        uint8_t j = sigma[i];
        
        std::cout << "(" << (uint32_t)i;
        while( j != i ) // mark all symbols in this cycle
        {
            std::cout << " " << (uint32_t)j;
            visited[j] = true;
            j = sigma[j];
        }
        std::cout << ")";

        // find the next unvisited cycle
        for( ++i; i <= p && visited[i]; ++i )
        {
        }
    }
}

size_t HashTuple :: operator ()( const Tuple &tuple ) const
{
    size_t hashvalue=0;
    size_t offset = 2;
    for( const auto& cit : tuple.rep )
    {
        hashvalue += offset*(cit.first + 8*cit.second);
        offset *= 16;
    }
    return hashvalue;
}


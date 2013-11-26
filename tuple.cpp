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
    #ifdef KAPPA_DEBUG_TUPEL
    PermutationType pt = permutation_type();
    stream << "Externe = " << pt.num_cycles - 1 << " Interne = " << pt.num_punctures << " p = " << p << ": ";
    #endif
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
    for( uint32_t i = 1; i < norm() - 1; i++ )
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
        auto a = boundary[q].first;
        auto b = boundary[q].second;
        
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

Tuple Tuple :: d_hor_test( uint8_t i ) const
{
    if( (i == 0) || (i >= this->norm()) )
    {
        return Tuple(this->norm());
    }
    Tuple boundary = *this;
    Tuple::Permutation sigma = long_cycle();     // sigma
    Tuple::Permutation sigma_inv = long_cycle_inv(); // sigma^{-1}

    // Compute the permutations sigma_i. The only symbols that change when we multiply with tau_i = (k,l) are
    // sigma_{i-1}^{-1}(k) and sigma_{i-1}^{-1}(l).
    for( uint32_t l = 1; l <= boundary.norm(); l++ )
    {
        // write tau_l = (a,b)
        auto& a = boundary.at(l).first;
        auto& b = boundary.at(l).second;
        
        // compute sigma_l = tau_l sigma_{l-1}
        std::swap( sigma[ sigma_inv[ a ] ], sigma[ sigma_inv[ b ] ] );

        // compute sigma_l^{-1}
        std::swap( sigma_inv[ a ], sigma_inv[ b ] );

        if( a == i )  // This row is to be deleted.
        {
            a = sigma[ a ]; // With Lemma 66 one sees that only the symbol i -> sigma(i) has to be changed.
            if( a == i )  // The row is degenerate.
            {
                return Tuple( boundary.norm() );
            }
        }
        if( a > i )   // all rows above i decrease
        {
            --a;
        }

        if( b == i ) // similar as above
        {
            b = sigma[ b ];
            if( b == i)
            {
                return Tuple( boundary.norm() );
            }
        }
        if( b > i )
        {
            --b;
        }

        if( a < b )    // The bigger symbol always is to the left of the other one.
        {
            std::swap(a,b);
        }
    }

    return boundary;
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

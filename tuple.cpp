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

Transposition& Tuple :: operator[](size_t n)
{
    return rep[n-1];
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
    uint32_t i = 0;
    for( std::vector< Transposition >::const_reverse_iterator it = tuple.rep.crbegin(); it != tuple.rep.crend(); ++it )
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
        ++i;
    }
    
    return stream;
}

PermutationType Tuple :: permutation_type()
{
    // Man kann ebenso die Anzahl der Zykel von (p p-1 ... 1) t_1 ... t_h zaehlen, da die t_i Transpositionen sind.
    PermutationType pt;    // Anzahl.first enstspricht den externen Punktierungen, Anzahl.second den internen.
    std::map< uint8_t, uint8_t > sigma_inv;

    // Initialisiere simga^{-1} als (p p-1 ... 1)
    for( uint8_t i = p; i > 1; i-- )
    {
        sigma_inv[i] = i-1;
    }
    sigma_inv[1] = p;

    // Multipliziere mit t_1 bis t_h und zaehle die internen Punktierungen
    for( uint32_t i = 0; i < rep.size(); i++ )
    {
        if(rep[i].first == rep[i].second)
        {
            pt.num_punctures += 1;
        }

        uint8_t k = rep[i].first;
        uint8_t l = rep[i].second;
        std::swap( sigma_inv[ k ], sigma_inv[ l ] );
    }

    // Zaehle die Zylel
    std::vector<bool> besucht( p+1, false );    // besucht[0] wird gar nicht verwendet.
    for( uint8_t i = 1; i <= p; ) // Wir verarbeiten alle Zykel nacheinander und markieren, welche Symbole verwendet werden.
    {
        // Verarbeite den naechsten Zykel
        pt.num_cycles += 1;
        besucht[i] = true;
        uint8_t j = sigma_inv[i];
        
        while( j != i ) // Markiere alle Symbole dieses Zykels.
        {
            besucht[j] = true;
            j = sigma_inv[j];
        }

        // Finde den naechsten unbesuchten Zykel
        for( ++i; i <= p && besucht[i]; ++i )
        {
        }
    }
    return pt;
}

bool Tuple :: monoton()
{
    // Ein Tupel ist genau dann monoto, wenn die Folge rep[i] monoton ist.
    for( uint32_t i = 0; i < rep.size() - 1; i++ )
    {
        if( rep[i+1].first < rep[i].first )
        {
            return false;
        }
    }

    return true;
}

bool Tuple :: f(uint32_t i)
{
    // Compare Mehner: page 49.
    // Attention: The i-th Transposition is stored in rep[i-1].
    // Wird durch Fallunterscheidung klar.
    if( i >= rep.size() ) // i muss kleiner als h sein.
    {
        std::cerr << "Error in 'bool Tuple :: f(uint32_t i)' -> i=" << i << " >= h=" << rep.size() << std::endl;
        return false;
    }
    else    // (ab)(cd)
    {
        uint8_t a = rep[i].first;
        uint8_t b = rep[i].second;
        uint8_t c = rep[i-1].first;
        uint8_t d = rep[i-1].second;
        
        if( a == c ) // (a )(a )
        {
            if( b == d ) // (ab)(ab) = id und (a*)(a*) = 0
            {
                rep.assign( rep.size(), Transposition(0,0) );
                return false;
            }
            else if( a == b ) // (a*)(ad) = (ad*) -> (d*)(ad)
            {
                rep[i] = Transposition(d,d);
                return true;
            }
            else if( a == d ) // (ab)(a*) = (ab*) -> (b*)(ab)
            {
                rep[i]   = Transposition(b,b);
                rep[i-1] = Transposition(a,b);
                return true;
            }
            else // (ab)(ad) = (adb) = (bd)(ab) bzw (db)(ab)
            {
                rep[i]   = Transposition(std::max(b,d), std::min(b,d));
                rep[i-1] = Transposition(a,b);
                return true;
            }
        }
        else if( a == b ) // (a*)(  )
        {
            if( a == d ) // (a*)(ca) = (ca*) -> (a*)(ca)
            {
                return true;
            }
            else // (a*)(c*) oder (a*)(cd) und beide kommutieren.
            {
                if( a < c )
                {
                    return true;
                }
                else
                {
                    std::swap( rep[i], rep[i-1]);
                    return true;
                }
            }
        }
        // Case: a > b
        else if( b == c ) // (ab)(b )
        {
            if( c == d ) // (ab)(b*) = (ab*) -> (b*)(ab)
            {
                std::swap( rep[i], rep[i-1]);
                return true;
            }
            else // (ab)(bd) = (abd) und (abd)' = (ad) und (abd)(ad) = (bd)
            {
                rep[i]   = Transposition(b,d);
                rep[i-1] = Transposition(a,d);
                return true;
            }
            // These are all cases since a > b >= d.
        }
        else if( a == d ) // (ab)(ca) = (acb) und c > a -> (ab)(ca)
        {
            return true;
        }
        else if( b == d ) // (ab)(cb) = (abc) -> a < c: (ab)(cb) oder a > c:(cb)(ac) 
        {
            if( a < c )
            {
                return true;
            }
            else    // (abc)' = (ac) und (abc)(ac) = (cb)
            {
                rep[i]   = Transposition(c,b);
                rep[i-1] = Transposition(a,c);
                return true;
            }
        }
        else // (ab)(c ) mit d nicht a oder b, also disjunkte Transpositionen.
        {
            if( a < c )
            {
                return true;
            }
            else
            {
                std::swap( rep[i], rep[i-1]);
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

    for( uint32_t j = q-1; j >= i; j-- ) // Da i > 0 erfuellt ist, endet die Schleife wirklich.
    {
        #ifdef KAPPA_DEBUG_TUPEL
        std::cerr << "    f_" << j << "( " << *this << " ) = ";
        #endif
        if( f(j) == false ) // Die Norm des Produkts faellt.
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

Tuple Tuple :: d_hor( uint8_t k )
{
    if(0 == k || k >= p)
    {
        return Tuple();
    }
    
    Tuple boundary = *this;
    
    // start with sigma_0.
    Tuple::Permutation sigma = long_cycle();
    Tuple::Permutation sigma_inv = long_cycle_inv();
    
    for(uint8_t q = 1; q <= boundary.rep.size(); ++q)
    {
        // Write tau_q = (a,b)
        auto a = boundary.rep[q].first;
        auto b = boundary.rep[q].second;
        
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
                    boundary.rep[q].first  = std::max(a,l);
                    boundary.rep[q].second = std::min(a,l);
                }
                else
                {
                    boundary.rep[q].first  = std::max(b,l);
                    boundary.rep[q].second = std::min(b,l);
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
    for(uint8_t q = 1; q <= boundary.rep.size(); ++q)
    {
        // Write tau_q = (a,b)
        auto a = boundary.rep[q].first;
        auto b = boundary.rep[q].second;
        
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

Tuple::Permutation Tuple::long_cycle()
{
    Tuple::Permutation sigma;
    for(uint8_t k = 1; k < p; ++k)
    {
        sigma[k] = k+1;
    }
    sigma[p] = 1;
    return sigma;
}

Tuple::Permutation Tuple::long_cycle_inv()
{
    Tuple::Permutation sigma;
    for(uint8_t k = 2; k <= p; ++k)
    {
        sigma[k] = k-1;
    }
    sigma[1] = p;
    return sigma;
}

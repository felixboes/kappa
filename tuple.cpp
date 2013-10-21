#include "tuple.hpp"



Tuple :: Tuple(size_t h) :
    p(0),
    rep( h, Transposition(0, 0) )
{
}

Transposition& Tuple :: operator[](size_t n)
{
    return rep[n];
}

int32_t Tuple :: norm() const
{
    if( *this == true )
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
    for( auto it = rep.rbegin(); it != rep.rend(); ++it )
    {
        stream << "(" << (uint32_t)it->first << ",";
        if( it->first != it->second )
            stream << (uint32_t)it->second << ")";
        else
            stream << "*)";
    }
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
    return anzahl;
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
    // Wird durch Fallunterscheidung klar.
    if( i >= rep.size() ) // i muss kleiner als h sein.
    {
        std::cerr << "f(i) mit i=" << i << "und h=" << rep.size() << std::endl;
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
            else // (ab)(ac) = (acb) und (acb)' = (ab) und (acb)(ab) = (bc) bzw. (cb)
            {
                rep[i]   = Transposition(a,b);
                rep[i-1] = Transposition(std::max(b,c), std::min(b,c));
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
}

bool Tuple :: phi( uint32_t q, uint32_t i, std::stringstream *s )
{
    if( i == 0 || i > q )
    {
        std::cerr << "phi^q_i mit q=" << q << "und i=" << i << std::endl;
        return false;
    }

    if( i == q ) // \Phi^q_q ist die Identitaet.
    {
        return true;
    }

    for( uint32_t j = q-1; j >= i; j-- ) // Da i > 0 erfuellt ist, endet die Schleife wirklich.
    {
        #ifdef KAPPA_DEBUG_TUPEL
        if(s != NULL)
        {
            *s << "    f_" << j << "( " << to_string() << " ) = ";
        }
        #endif
        if( f(j) == false ) // Die Norm des Produkts faellt.
        {
            #ifdef KAPPA_DEBUG_TUPEL
            if( s != NULL )
            {
                s->str( std::string() );
            }
            #endif
            return false;
        }
        #ifdef KAPPA_DEBUG_TUPEL
        else if ( s != NULL)
        {
            *s << to_string() << std::endl;
        }
        #endif
    }
    return true;
}

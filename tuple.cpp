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

bool Tuple :: phi( uint32_t q, uint32_t i )
{
    if( i == 0 || i > q )
    {
        std::cerr << "phi^q_i mit q=" << q << "und i=" << i << std::endl;
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

Tuple Tuple :: del2(uint32_t j)
{
    Tuple rand = *this;
//    std::map< uint8_t, uint8_t > sigma;     // sigma
//    std::map< uint8_t, uint8_t > sigma_inv; // sigma^{-1}

//    // initialisiere simga als (1 2 ... p)
//    sigma[1] = 2;
//    sigma_inv[1] = p;
//    for( uint8_t l = 2; l < p; l++ )
//    {
//        sigma[l] = l+1;
//        sigma_inv[l] = l-1;
//    }
//    sigma[p] = 1;
//    sigma_inv[p] = p-1;

//    // Berechne nun die sigma_i. Die beiden einzigen Symbole, die sich aendern, wenn wir mit tau_i = (k,l) mulitplizieren, sind
//    // sigma_{i-1}^{-1}(k) und sigma_{i-1}^{-1}(l).
//    for( uint32_t l = 1; l <= rep.size(); l++ )
//    {
//        if( rand.rep[l-1].first == j )  // Diese Zeile soll geloescht werden.
//        {
//            // Mit Lemma 66 sieht man ein, dass nur ein Symbol geaendert werden muss und zwar j -> sigma_{i-1}(j)
//            // Die Zelle ist genau dann degeneriert, wenn j = sigma_{i-1}(j) oder tau_i = (j, sigma_{i-1}(j)) für den unpunktierten Fall
//            // Für den punktierten Fall bekommen muss zusätzlich geprüft werden, ob tau_i = (j,*) und ?

//            rand.rep[l-1].first = sigma[ rand.rep[l-1].first ];
//            if( rand.rep[l-1].first == j || rand.rep[l-1].first == rand.rep[l-1].second )
//            {
//                rand = Tupel( rep.size() );
//                return rand;
//            }
//        }
//        if( rand.rep[l-1].second == j ) // Analog zu obigem Abschnitt.
//        {
//            rand.rep[l-1].second = sigma[ rand.rep[l-1].second ];
//            if( rand.rep[l-1].second == j || rand.rep[l-1].first == rand.rep[l-1].second)
//            {
//                rand = Tupel(rep.size());
//                return rand;
//            }
//        }
//        if( rand.rep[l-1].first > j )   // Alle Zeilen ueber i muessen nachrutschen.
//        {
//            rand.rep[l-1].first -= 1;
//        }
//        if( rand.rep[l-1].second > j )
//        {
//            rand.rep[l-1].second -= 1;
//        }

//        if( rand.rep[l-1].first < rand.rep[l-1].second )    // Das groessere Symbol soll immer links stehen.
//        {
//            uint8_t tmp = rand.rep[l-1].first;
//            rand.rep[l-1].first = rand.rep[l-1].second;
//            rand.rep[l-1].second = tmp;
//        }

//        // Berechne sigma_i = tau_i sigma_{i-1}
//        sigma[ sigma_inv[ rep[l-1].first ] ] = rep[l-1].second;
//        sigma[ sigma_inv[ rep[l-1].second ] ] = rep[l-1].first;

//        // Berechne sigma_i^{-1}
//        uint8_t tmp = sigma_inv[ rep[l-1].first ];
//        sigma_inv[ rep[l-1].first ] = sigma_inv[ rep[l-1].second ];
//        sigma_inv[ rep[l-1].second ] = tmp;

//    }

    std::cerr << " TODO ^^ " << std::endl;
    return rand;
}

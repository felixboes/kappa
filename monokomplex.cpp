#include "monokomplex.hpp"

/*
 *
 *   MonoModul
 *
 */

void MonoModul :: show()
{
    for( std::vector< Tuple >::iterator it = omega.begin(); it != omega.end(); ++it )
    {
        it->show();
        std::cout << " ";
    }
}

uint64_t MonoModul :: size()
{
    return omega.size();
}

void MonoModul :: ergaenze (Tuple& erzeuger)
{
    omega.push_back(erzeuger);
}

/*
 *
 *   MonoKomplex
 *
 */

MonoKomplex :: MonoKomplex(uint32_t _g, uint32_t _m, bool _externe) : g(_g), m(_m), h(2*_g + _m), externe(_externe)
{
    Tuple anf(h);

    if(externe)
    {
        anf[0] = std::pair< uint8_t, uint8_t >(2, 1);
        gen_modules_extern(1, 2, anf);  // Wir beginnen mit ...(2 1)
    }
    else
    {
        anf[0] = std::pair< uint8_t, uint8_t >(1, 1);
        gen_modulen_intern(1, 1, anf, 1); // Wir beginnen mit einer internen Punktierung ...(1 1)

        anf[0] = std::pair< uint8_t, uint8_t >(2, 1);
        gen_modulen_intern(1, 2, anf, 0); // Wir beginnen mit ...(2 1)
    }

}

void MonoKomplex :: show()
{

    for( std::map< uint32_t, MonoModul >::iterator it = modul.begin(); it != modul.end(); ++it )
    {
        std::cout << "Modul " << h + it->first << std::endl;
        it->second.show();
        std::cout << std::endl;
    }
}

void MonoKomplex :: show_size()
{
    std::cout << "Ich liste die Groesse der Moduln der obersten Zeile auf. Dabei ist der horizontale Grad angegeben." << std::endl;
    uint64_t gesamt = 0;
    for( std::map< uint32_t, MonoModul >::iterator it = modul.begin(); it != modul.end(); ++it )
    {
        gesamt += it->second.size();
        std::cout << "Grad p=" << it->first << ": " << it->second.size() << std::endl;
    }
    std::cout << "Gesamt: " << gesamt << std::endl;
}

void MonoKomplex :: gen_modules_extern(uint32_t l, uint32_t p, Tuple& tuple)
{
    /* Up to now we have determined all monotonic tuples of l-1 transpositions containing the 
       symbols 1, ..., p, each at least once. We now add an l-th transposition and continue
       recursively.*/
    if(l < h) // There are h-l transpositions left to be determined.
    {
	/* From an (l-1)-tuple containing p symbols we can build up an l-tuple 
        with p, p+1 or p+2 symbols. */

        /* p -> p
           In this case we use the same number of symbols. Since we only
           enumerate monotonic tuples, the height of the l-th transposition needs to be p.
           We try out all possibilities for the second symbol in the l-th transposition. */
        tuple.p = p;
        for(uint32_t i = p-1; i > 0; i--)
        {
            tuple[l] = std::pair< uint8_t, uint8_t >(p, i); 
            gen_modules_extern(l+1, p, tuple);
        }

        /* p -> p+1
           In this case we use p+1 symbols instead of p. Two cases occur. */
        tuple.p = p+1;
        /* Case 1: The new row in the Schlitzgebiet is inserted at the top, i.e. 
                   the transpositions 1, ..., l-1 remain the same and we only insert the symbol 
                   p+1 in the l-th transposition, together with any symbol of 1, ..., p. */
        for(uint32_t i = p; i > 0; i--)
        {
            tuple[l] = std::pair< uint8_t, uint8_t >(p+1, i);
            gen_modules_extern(l+1, p+1, tuple);
        }

        /* Case 2: The new row in the Schlitzgebiet is not inserted at the top but at a position
                   i = 1, ..., p. Thus all indices i+1, ..., p are shifted up by one. 
                   TODO: Check the indices!*/
        for( uint32_t i = p; i > 0; i-- )
        {
            Tuple tmp = tuple;
            for( uint32_t j = l; j > 0; j-- )
            {
                if( tmp[j-1].first >= i )
                {
                    tmp[j-1].first++;
                    if(tmp[j-1].second >= i )
                    {
                        tmp[j-1].second++;
                    }
                }
            }

            tmp[l] = std::pair< uint8_t, uint8_t >(p+1, i);
            gen_modules_extern(l+1, p+1, tmp);
        }

        // p -> p+2
        tuple.p = p+2;
        for( uint32_t i = p+1; i > 0; i-- )
        {
            Tuple tmp = tuple;
            for( uint32_t j = l; j > 0; j-- )
            {
                if( tmp[j-1].first >= i )
                {
                    tmp[j-1].first++;
                    if(tmp[j-1].second >= i )
                    {
                        tmp[j-1].second++;
                    }
                }
            }

            tmp[l] = std::pair< uint8_t, uint8_t >(p+2, i);
            gen_modules_extern(l+1, p+2, tmp);
        }
    }
    else
    {
        if(tuple.zykelzahl().first == m+1)
        {
            modul[p].ergaenze( tuple );
        }
    }
}

void MonoKomplex :: gen_modulen_intern(uint32_t l, uint32_t p, Tuple& tuple, uint32_t intpkt)
{
    if( intpkt <= m && h - l >= m - intpkt ) // Die Anzahl unverwendeter Transpositionen muss mindestens so gross sein, wie die Anzahl noch zu vergebenden gefaerbten Fixpunkte.
    {
        if( l < h )
        {
            // p -> p
            tuple.p = p;

            if( intpkt < m ) // Es sind interne Punktiereungen uebrig
            {
                tuple[l] = std::pair< uint8_t, uint8_t >(p, p);
                gen_modulen_intern(l+1, p, tuple, intpkt+1);
            }
            for(uint32_t i = p-1; i > 0; i--)
            {
                tuple[l] = std::pair< uint8_t, uint8_t >(p, i);
                gen_modulen_intern(l+1, p, tuple, intpkt);
            }

            // p -> p+1
            // die neue Zeile ist oben.
            tuple.p = p+1;
            if( intpkt < m ) // Es sind interne Punktiereungen uebrig
            {
                tuple[l] = std::pair< uint8_t, uint8_t >(p+1, p+1);
                gen_modulen_intern(l+1, p+1, tuple, intpkt+1);
            }
            for(uint32_t i = p; i > 0; i--)
            {
                tuple[l] = std::pair< uint8_t, uint8_t >(p+1, i);
                gen_modulen_intern(l+1, p+1, tuple, intpkt);
            }

            // p -> p+1
            // die neue Zeile ist nicht oben
            for( uint32_t i = p; i > 0; i-- )
            {
                Tuple tmp = tuple;
                for( uint32_t j = l; j > 0; j-- )
                {
                    if( tmp[j-1].first >= i )
                    {
                        tmp[j-1].first++;
                        if(tmp[j-1].second >= i )
                        {
                            tmp[j-1].second++;
                        }
                    }
                }

                tmp[l] = std::pair< uint8_t, uint8_t >(p+1, i);
                gen_modulen_intern(l+1, p+1, tmp, intpkt);
            }

            // p -> p+2
            tuple.p = p+2;
            for( uint32_t i = p+1; i > 0; i-- )
            {
                Tuple tmp = tuple;
                for( uint32_t j = l; j > 0; j-- )
                {
                    if( tmp[j-1].first >= i )
                    {
                        tmp[j-1].first++;
                        if(tmp[j-1].second >= i )
                        {
                            tmp[j-1].second++;
                        }
                    }
                }

                tmp[l] = std::pair< uint8_t, uint8_t >(p+2, i);
                gen_modulen_intern(l+1, p+2, tmp, intpkt);
            }
        }
        else
        {
            if( intpkt == m && tuple.zykelzahl().first == 1 )
            {
                modul[p].ergaenze( tuple );
            }
        }
    }
}

void MonoKomplex :: del2_kappa(uint32_t p)
{
    /**
        Statt die Funktion rekusiv aufzurufen, zaehlen wir die Indexfolgen ab.
        Dazu sieht man ein, dass die Abbildung
        \f[
            \{ 1, \ldots, n! \} \to \{ (s_{n-1}, \ldots, s_0) \mid 0 \le s_q < q \} \hspace{15pt} k \mapsto
            \left( \left\lfloor \frac{k}{q!} \right\rfloor \mathrm{mod}\,q \right)_{q}
        \f]
        eine Bijektion ist.
    **/

    uint64_t mono = 0;

    std::cout << "Ich berechne nun del o kappa_" << p << " ... ";
    #ifdef KAPPA_DEBUG_MONOKOMPLEX
         std::cout << std::endl;
    #endif

    // Fuer jedes Tuple berechnen wir alle in Kappa vorkommenden Basiselemente.
    for( std::vector< Tuple >::iterator it = modul[p].omega.begin(); it != modul[p].omega.end(); ++it )
    {
        #ifdef KAPPA_DEBUG_MONOKOMPLEX
        std::cout << std::setfill('-') << std::setw(40) << "-" << std::setfill(' ') << std::endl
                  << "Betrachte Tuple: " << it->to_string() << std::endl
                  << std::setfill('-') << std::setw(40) << "-" << std::setfill(' ') << std::endl;
        #endif

        // Um Variablen mit openmp fuer jeden Thread einzeln zu hinterlegen, muessen sie vor der Schleife definiert werden.
        Tuple neuer;
        Tuple rand;
        bool valid;
        uint32_t k;
        uint32_t i;
        uint32_t q;
        uint32_t s_i;
        #ifdef KAPPA_PARA
        #pragma omp parallel for private(neuer,rand,valid,k,i,q,s_i) shared(mono)
        #endif
        for( k = 0; k < factorial(h); k++ )        // Abzaehlung fuer die Tuple (s_n, ..., s_1) mit 1 <= s_i <= i
        {
            neuer = *it;
            valid = true;

            #ifdef KAPPA_DEBUG_MONOKOMPLEX
            std::stringstream f_folge;
            std::string indexfolge = "1)";
            #endif

            for( q = 2; q <= h; q++ )
            {
                s_i = 1 + ( ( k / factorial(q-1)) % q );    // Abzaehlung fuer die Tuple (s_h, ..., s_1) mit 1 <= s_i <= i

                #ifdef KAPPA_DEBUG_MONOKOMPLEX
                std::stringstream tmp;
                tmp << s_i;
                indexfolge = tmp.str() + "," + indexfolge;

                if( neuer.phi(q, s_i, &f_folge) == false )
                {
                    valid = false;
                    f_folge.str( std::string() );
                    break;
                }
                #else
                if( neuer.phi(q, s_i) == false )
                {
                    valid = false;
                    break;
                }
                #endif
            }
            if( valid )   // Berechne alle horizontalen Raender
            {
                #ifdef KAPPA_DEBUG_MONOKOMPLEX
                indexfolge = "(" + indexfolge;
                std::cout << "Phi_" << indexfolge << "( " << it->to_string() << " ) = " << neuer.to_string() << std::endl;
                if( f_folge.str().length() > 1 )
                {
                    std::cout << "Dabei ist:" << std::endl
                              << f_folge.str();
                }
                f_folge.str( std::string() );
                bool rand_vorhanden = false;
                #endif
                for( i = 1; i < p; i++ )
                {
                    rand = neuer.del2(i);
                    if( rand.() )
                    {
                        #ifdef KAPPA_DEBUG_MONOKOMPLEX
                        rand_vorhanden = true;
                        f_folge << "    Der " << i << "-te Rand " << rand.to_string();
                        #endif
                        if( rand.monoton() == true ) // Traegt zum Differential bei.
                        {
                            #ifdef KAPPA_PARA
                            // Wenn wir mit openmp parallelisieren, muss mono++ atomar sein.
                            #pragma omp atomic
                            #endif
                            mono++;
                            #ifdef KAPPA_DEBUG_MONOKOMPLEX
                            f_folge << " ist monoton." << std::endl;
                            #endif
                        }
                        #ifdef KAPPA_DEBUG_MONOKOMPLEX
                        else
                        {
                            f_folge << " ist nicht monoton." << std::endl;
                        }
                        #endif
                    }
                }
                #ifdef KAPPA_DEBUG_MONOKOMPLEX
                if( rand_vorhanden == true )
                {
                    std::cout << "Die horizontalen Raender sind:" << std::endl
                              << f_folge.str() << std::endl;
                }
                else
                {
                    std::cout << "Alle horizontalen Raender sind degeneriert." << std::endl << std::endl;
                }
                #endif
            }
        }
    }

    #ifdef KAPPA_DEBUG_MONOKOMPLEX
    std::cout << "Die Anzahl der monotonen Paare (t,t',i), wobei t' = del_i o kappa (t) ist, betraegt genau " << mono << "." << std::endl;
    #else
    std::cout << mono << std::endl;
    #endif
}

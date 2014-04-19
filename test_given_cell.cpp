#include <iostream>
#include <string>

#include "factorial.hpp"
#include "tuple.hpp"

//          Zelle aus Text          //
//          Zelle aus Text          //
//          Zelle aus Text          //
struct Token
{
    enum Typ {
        kl_auf,
        zahl,
        komma,
        kl_zu,
        bar,
        ende,
        error
    };

    Token( Typ _t = ende, uint32_t _z = 0) : t(_t), z(_z) {}

    Typ t;
    uint8_t z;

    operator bool() { return ( t != ende && t != error ); }
};

class Parser
{
public:
    Parser (char* s) : str(s) {}

    // Lies nächstes Token.
    Token operator()()
    {
        Token t;
        if( str == nullptr )
            return Token( Token::error, 0 );
    
        if( *str == '(' )
            t = Token( Token::kl_auf, 0 );
        else if ( *str == ')')
            t = Token( Token::kl_zu, 0 );
        else if ( *str == ',')
            t = Token( Token::komma, 0 );
        else if ( *str == '|')
            t = Token( Token::bar, 0 );
        else if ( *str >= '0' && *str <= '9' )
        {
            // Lese Zahl
            t = Token( Token::zahl, atoi(str) );
            // Gehe an letzte Stelle der Zahl
            while( *str >= '0' && *str <= '9' )
            {
                str++;
            }
            return t;
        }
        else if ( *str == '\0' )
            t = Token( Token::ende, 0 );
        else // Versuche nächstes Zeichen zu lesen
        {
            str++;
            t = this->operator()();
        }
        str++;
        return t;
    }

private:
    char* str;
};

Tuple parse( int32_t argc, char** argv )
{
    if(argc < 4)
    {
        return Tuple();
    }

    uint32_t sym = atoi(argv[1]);
    uint32_t h = atoi(argv[2]);

    Tuple zelle = Tuple(sym,h);

    Token::Typ state = ( Token::kl_auf );
    uint32_t k = h;

    // Iteriere duch die ggf. gestückelte Zellenschreibweise
    for(int32_t i = 3; i <= argc; ++i )
    {
        Parser p = Parser(argv[i]);
        Token t;
        while( t = p() )
        {
            //std::cout << "State: " << state << " Token: " << t.t << " k: " << k << " Zahl: " << (int32_t)t.z << std::endl;
            switch (state)
            {

            case Token::kl_auf:
                if( t.t == Token::kl_auf )
                {} // Weiteres '(' ignorieren.
                else if( t.t == Token::zahl )
                {
                    zelle[k].first = t.z;
                    state = Token::zahl;
                }
                else
                {
                    return Tuple(); // Zelle konnte nicht eingelesen werden.
                }
                break;

            case Token::zahl:
                if( t.t == Token::komma )
                {
                    state = Token::komma;
                }
                else if ( t.t == Token::kl_zu )
                {
                    state = Token::kl_zu;
                }
                else
                {
                    return Tuple();
                }
                break;

            case Token::komma:
                if( t.t == Token::zahl )
                {
                    zelle[k].second = t.z;
                    state = Token::zahl;
                }
                else
                {
                    return Tuple();
                }
                break;

            case Token::kl_zu:
                if( t.t == Token::kl_zu )
                {} // Weiteres ')' ignorieren.
                else if( t.t == Token::bar ) // Nächste Transposition
                {
                    k--;
                    state = Token::bar;
                }
                else
                {
                    return Tuple();
                }
                break;

            case Token::bar:
                if( t.t == Token::kl_auf )
                {
                    state = Token::kl_auf;
                }
                else
                {
                    return Tuple();
                }
                break;

            case Token::error:
                break; // Konnte Zeichen nicht lesen. Weitermachen

            case Token::ende:
                break; // String ist am Ende.

            default:
                return Tuple();
                break;
            }
        }
    }

    return zelle;
}
//          Zelle aus Text          //
//          Zelle aus Text          //
//          Zelle aus Text          //

//          Brechstange             //
//          Brechstange             //
//          Brechstange             //
void S_(size_t i, Tuple& zelle )
{
    auto k = zelle[i].first;
    for( auto j = i+1; j <= zelle.norm(); ++j )
    {
        zelle[j].first++;
        if( zelle[j].second >= k )
        {
            zelle[j].second++;
        }
    }
    zelle.p++;
}
//          Brechstange             //
//          Brechstange             //
//          Brechstange             //

std::string s_q_folge(uint32_t k, uint32_t h)
{
    std::string str = "1)";
    for( uint32_t q = 2; q <= h; q++ )
    {
        str = std::to_string( 1 + ( ( k / factorial(q-1)) % q ) ) + "," + str;
    }
    
    return "(" + str;
}

int main( int argc, char** argv )
{
    Tuple zelle = parse(argc, argv);
    S_(1,zelle);
    std::cout << zelle << std::endl;
    
    uint32_t h = zelle.norm();
    uint32_t p = zelle.p;
    
    for( uint32_t k = 0; k < factorial(h); k++ )
    {
        Tuple test_zelle = zelle;
        uint32_t s_q;
        bool norm_preserved = true;
        
        // Calculate phi_{(s_h, ..., s_1)}( Sigma )
        for( uint32_t q = 1; q <= h; q++ )
        {
            s_q = 1 + ( ( k / factorial(q-1)) % q );   
            if( test_zelle.phi(q, s_q) == false )
            {
                norm_preserved = false;
                break;
            }
        }
        
        // Calculate its faces
        if( norm_preserved )   // Compute all horizontal faces.
        {
            Tuple boundary;
            std::cout << "K_" << s_q_folge(k,h) << " -> " << test_zelle << std::endl;
            for( uint32_t i = 1; i < p; i++ )
            {
                if( (boundary = test_zelle.d_hor(i)) && boundary.monotone() )
                {
                    std::cout << "           d_" << i << " -> " << boundary << std::endl;
                }
            }
        }
        std::cout << std::endl;
    }
    
    return 0;
}

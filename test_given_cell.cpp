#include <iostream>
#include <string>

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




int main( int argc, char** argv )
{
    std::cout << parse(argc, argv) << std::endl;
    return 0;
}

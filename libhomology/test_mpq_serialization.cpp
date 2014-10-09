#include <cstdio>
#include <iostream>
#include <string>

#include <gmpxx.h>

int main( int argc, char** argv )
{
    if( argc == 1 || argv[1] == "c" )
    {
        std::string bz2_command_string = std::string("bzip2 -c -z -7 > test.bz2");
        const char* bz2_command = bz2_command_string.c_str();
        FILE* bz2_pipe;
        
        if( (bz2_pipe = popen(bz2_command, "w")) == nullptr )
        {
            std::cout << "Fehler beim Öffnen." << std::endl;
            return 1;
        }
        
        mpq_class t=200;
        t /= 7;
        
        std::cout << t << std::endl;
        
        mpz_out_raw( bz2_pipe, t.get_num().get_mpz_t() );
        mpz_out_raw( bz2_pipe, t.get_den().get_mpz_t() );
        
        pclose(bz2_pipe);
        return 0;
    }
    else
    {
        std::string bz2_command_string = std::string("bzip2 -c -d < test.bz2");
        const char* bz2_command = bz2_command_string.c_str();
        FILE* bz2_pipe;
        
        if( (bz2_pipe = popen(bz2_command, "r")) == nullptr )
        {
            std::cout << "Fehler beim Öffnen." << std::endl;
            return 1;
        }
        
        mpz_t in;
        mpz_init(in);
        while( mpz_inp_raw( in, bz2_pipe ) != 0 )
        {
            std::cout << in << std::endl;
        }
        
        pclose(bz2_pipe);
        return 0;
    }
}

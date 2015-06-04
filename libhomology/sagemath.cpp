#include "sagemath.hpp"

template<>
std::string sage_coeff<Q>()
{
    return "QQ";
}

SagemathInterface::SagemathInterface()
{
    // Open sage.
    std::string sagemath_command_string = "cat > ./sagecmds";
    std::cout << "calling " << sagemath_command_string << std::endl;
    const char* sagemath_command = sagemath_command_string.c_str();

    if( (sagemath_pipe = popen( sagemath_command, "w")) == nullptr )
    {
        std::cout << "Error: Could not open sagemath interface." << std::endl;
        pclose(sagemath_pipe);
        return;
    }

    fprintf( sagemath_pipe, "f = file('output.txt','w')\n" );
    fprintf( sagemath_pipe, "A=0\n" );
}

SagemathInterface::~SagemathInterface()
{
    fprintf( sagemath_pipe, "f.close()\n" );
    pclose(sagemath_pipe);
}

void SagemathInterface::test() const
{
    if( sagemath_pipe == nullptr )
    {
        std::cout << "pipe closed" << std::endl;
        return;
    }
    fprintf( sagemath_pipe, "A = matrix(QQ, 10, sparse=True)\n" );
    fprintf( sagemath_pipe, "A.str()\n" );
}

void SagemathInterface::compute_rank()
{
    if( sagemath_pipe == nullptr )
    {
        std::cout << "pipe closed" << std::endl;
        return;
    }
    fprintf( sagemath_pipe, "A.str()\n" );
    fprintf( sagemath_pipe, "A.rank()\n" );
    fprintf( sagemath_pipe, "rk = A.rank()\n" );
    fprintf( sagemath_pipe, "f.write(str(rk)+'\\n')\n" );
}

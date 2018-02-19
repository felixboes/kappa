#include <chrono>
#include <iostream>
#include <thread>

#include <boost/numeric/ublas/io.hpp>

#include "diagonalizer_zm.hpp"
#include "matrix_zm.hpp"

int main(int argc, char** argv)
{
    if( argc <= 1 )
    {
        std::cout << "Modulus nicht eingegeben. Beende." << std::endl;
        return 1;
    }
    
    Zm::set_modulus(atoi(argv[1]), 1);
    
    std::cout.setf(std::ios::unitbuf);
    
    std::cout << "Anzahl Zeilen eingeben." << std::endl;
    uint32_t zeilen;
    std::cin >> zeilen;
    
    std::cout << "Anzahl Spalten eingeben." << std::endl;
    uint32_t spalten;
    std::cin >> spalten;
    
    MatrixZm Test(zeilen,spalten);
    std::cout << Test.data().size() << std::endl;
    
    Test.resize(0,0);
    std::cout << Test.data().size() << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(5000));
    
    return 0;
}

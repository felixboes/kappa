#include <boost/numeric/ublas/io.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <chrono>
#include <cinttypes>
#include <functional>
#include <iostream>
#include <limits>

#include "chain_complex.hpp"
#include "diagonalizer_zm.hpp"
#include "homology_field.hpp"
#include "matrix_zm.hpp"

// This is a deterministic random generator.
// The initial seed is the system time.
boost::random::mt19937 gen_zm(std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()));

bool create_random_matrix_zm(uint32_t rows, uint32_t cols, uint32_t rank)
{
    MatrixZm matrix(rows, cols);
    MatrixZm row_ops = MatrixZmIdentity(rows);
    MatrixZm col_ops = MatrixZmIdentity(cols);
    
    if( rank > std::min(rows, cols) )
    {
        std::cerr << "rank > min(rows, cols)" << std::endl;
        return false;
    }
    
    for( uint32_t i = 0; i < rank; ++i )
    {
        matrix(i,i) = 1;
    }
    
    boost::random::uniform_int_distribution<> number(-100, 100);
    // row_ops = lower left diag.
    for( uint32_t i = 1; i < rows; ++i )
    {
        for( uint32_t j = 0; j <= i-1; ++j )
        {
            row_ops(i,j) = Zm( number(gen_zm) );
        }
    }
    
    // col_ops = upper right sdiag.
    for( uint32_t i = 0; i < cols-1; ++i )
    {
        for( uint32_t j = i+1; j < cols ; ++j )
        {
            col_ops(i,j) = Zm( number(gen_zm) );
        }
    }
    
    matrix = boost::numeric::ublas::prod( row_ops, matrix );
    matrix = boost::numeric::ublas::prod( matrix,  col_ops );
    DiagonalizerZm diagonalizer;
    
    return rank == diagonalizer.diag_field(matrix);
}

uint32_t test_rank_zm( uint32_t num_rounds, uint32_t max_num_rows, uint32_t max_num_cols, uint32_t max_num_rank = std::numeric_limits<uint32_t>::max() )
{
    uint32_t number_of_errors = 0;
    boost::random::uniform_int_distribution<> rnd_rows(5, max_num_rows);
    boost::random::uniform_int_distribution<> rnd_cols(5, max_num_cols);
    for( uint32_t round = 0; round < num_rounds; ++round )
    {
        uint32_t rows = rnd_rows(gen_zm);
        uint32_t cols = rnd_cols(gen_zm);
        boost::random::uniform_int_distribution<> rnd_rank( 5, std::min( std::min ( rows, cols ), max_num_rank ) );
        uint32_t rank = rnd_rank(gen_zm);
        
        std::cout << "Round: " << round << "\r";
        if( !create_random_matrix_zm(rows, cols, rank) )
        {
            std::cerr << "Error in Round: " << round << std::endl;
            number_of_errors++;
        }
    }
    return number_of_errors;
}

void test_some_chain_complex_zm()
{
    typedef ChainComplex< Zm, MatrixZm, DiagonalizerZm, HomologyField > ChainComplexZm;
    ChainComplexZm cc;
    
    MatrixZm M(4,4);
    M(0,0) = 1;
    M(0,1) = 1;
    M(0,2) = 2;
    M(0,3) = 3;

    M(1,0) = 0;
    M(1,1) = 0;
    M(1,2) = 0;
    M(1,3) = 0;

    M(2,0) = 1;
    M(2,1) = 0;
    M(2,2) = 5;
    M(2,3) = 7;

    M(3,0) = 0;
    M(3,1) = 0;
    M(3,2) = 0;
    M(3,3) = 0;

    MatrixZm N(4,2);
    N(0,0) = -15;
    N(0,1) =  -7;
    
    N(1,0) =   9;
    N(1,1) =   4;
    
    N(2,0) =   3;
    N(2,1) =   0;
    
    N(3,0) =   0;
    N(3,1) =   1;

    cc[1] = M;
    cc[0] = N;
    
    std::cout << "M:    " << cc[1] << std::endl;
    std::cout << "N:    " << cc[0] << std::endl;
    std::cout << "Prod: " << prod( cc[1], cc[0] ) << std::endl;
    
    HomologyField ho = cc.homology(0);
    std::cout << ho << std::endl;
}

int main( int argc, char ** argv )
{
    std::cout.setf(std::ios::unitbuf);
    
    if( argc <= 1 )
    {
        std::cerr << "Error: #Parameter = 0" << std::endl;
        return 1;
    }
    Zm::set_modulus(atoi(argv[1]),1);
    
    /*
    uint32_t errors = 0;
    if( argc == 5 )
    {
        errors = test_rank_zm( atoi(argv[2]), atoi(argv[3]), atoi(argv[4]) );
    }
    else if( argc >= 6 )
    {
        errors = test_rank_zm( atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]) );
    }
    else
    {
        std::cerr << "Error: #Parameter < 4" << std::endl;
        return 1;
    }
    std::cout << "Total number of errors: " << errors << std::endl;
    */
    test_some_chain_complex_zm();
    return 0;
}

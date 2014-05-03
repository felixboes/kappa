#include <boost/numeric/ublas/io.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <chrono>
#include <cinttypes>
#include <functional>
#include <iostream>

#include "chain_complex.hpp"
#include "diagonalizer_q.hpp"
#include "homology_field.hpp"
#include "matrix_q.hpp"


// This is a deterministic random generator.
// The initial seed is the system time.
boost::random::mt19937 gen_q(std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()));

/**
 * Creates a random matrix of size rows x cols of the given rank and tests where our
 * diagonalizer computes the rank of that matrix correctly.
 * \return true iff The rank of the computed matrix coincides with the rank of that matrix
 * computed by the diagonalizer.
 * Announces an error if there is no matrix of that size and rank.
 */
bool create_random_matrix_q(uint32_t rows, uint32_t cols, uint32_t rank)
{
    // To create a random matrix of the given rank, we start with a (non-random) matrix of
    // this rank. Then we multiply this matrix with invertible random matrices from the left to
    // obtain the desired random matrix.

    MatrixQ matrix(rows, cols);
    MatrixQ row_ops = MatrixQIdentity(rows);
    MatrixQ col_ops = MatrixQIdentity(cols);
    
    if( rank > std::min(rows, cols) )
    {
        std::cerr << "rank > min(rows, cols)" << std::endl;
        return false;
    }

    // Initialize matrix to be given by 1s on the first rank positions of the diagonal
    // and zeroes elsewhere.
    for( uint32_t i = 0; i < rank; ++i )
    {
        matrix(i,i) = 1;
    }
    
    // ranges for the random numbers used for the numerators and denominators
    boost::random::uniform_int_distribution<> numerator(-100, 100);
    boost::random::uniform_int_distribution<> denominator(1, 100);

    // Define the matrix row_ops as a random lower triangular matrix.
    for( uint32_t i = 1; i < rows; ++i )
    {
        for( uint32_t j = 0; j <= i-1; ++j )
        {
            row_ops(i,j) = Q( numerator(gen_q), denominator(gen_q) );
        }
    }
    
    // Define the matrix col_ops as a random upper triangular matrix.
    for( uint32_t i = 0; i < cols-1; ++i )
    {
        for( uint32_t j = i+1; j < cols ; ++j )
        {
            col_ops(i,j) = Q( numerator(gen_q), denominator(gen_q) );
        }
    }
    
    // Multiply the matrices to obtain a matrix of the given rank.
    matrix = boost::numeric::ublas::prod( row_ops, matrix );
    matrix = boost::numeric::ublas::prod( matrix,  col_ops );

    // Use the Diagonalizer to compute the rank of the random matrix.
    DiagonalizerQ diagonalizer;
    diagonalizer(matrix);
    
    return rank == diagonalizer.rank();
}

/**
 * Tests whether the diagonalizer computes the rank of matrices correctly via testing
 * it on various random matrices of random sizes.
 * \return number of matrices where the diagonalizer fails.
 * \param num_rounds number of tests
 * \param max_num_rows maximum number of rows of the test matrices
 * \param max_num_cols maximum number of columns of the test matrices
 * \max_rank maximum rank of the test matrices
 */
uint32_t test_rank_q( uint32_t num_rounds, uint32_t max_num_rows, uint32_t max_num_cols, uint32_t max_rank = std::numeric_limits<uint32_t>::max() )
{
    // variable used to count the failures of the diagonalizer
    uint32_t number_of_errors = 0;

    boost::random::uniform_int_distribution<> rnd_rows(5, max_num_rows);
    boost::random::uniform_int_distribution<> rnd_cols(5, max_num_cols);
    for( uint32_t round = 0; round < num_rounds; ++round )
    {
        uint32_t rows = rnd_rows(gen_q);
        uint32_t cols = rnd_cols(gen_q);
        boost::random::uniform_int_distribution<> rnd_rank( 5, std::min( std::min ( rows, cols ), max_rank ) );
        uint32_t rank = rnd_rank(gen_q);
        
        std::cout << "Round: " << round << "\r";
        if( !create_random_matrix_q(rows, cols, rank) )
        {
            std::cerr << "Error in Round: " << round << std::endl;
            number_of_errors++;
        }
    }
    return number_of_errors;
}

/**
 * Tests algorithms on chain complexes.
 */
void test_some_chain_complex_q()
{
    typedef ChainComplex< Q, MatrixQ, DiagonalizerQ, HomologyField > ChainComplexQ;
    ChainComplexQ cc;
    
    MatrixQ M(4,4);
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

    MatrixQ N(4,2);
    N(0,0) = -15;
    N(0,1) =  -7;
    
    N(1,0) =   9;
    N(1,1) =   4;
    
    N(2,0) =   3;
    N(2,1) =   0;
    
    N(3,0) =   0;
    N(3,1) =   1;

    cc[2] = N;
    cc[1] = M;
    cc[0] = MatrixQ(0,4);
    
    cc[10] = MatrixQ(0,1);
    cc[11] = MatrixQ(1,1);
    cc(11,0,0) = 0;
    cc[12] = MatrixQ(1,1);
    cc(12,0,0) = 2;
    cc[13] = MatrixQ(0,1);
    cc[14] = MatrixQ(1,1);
    cc(14,0,0) = 2;
    cc[15] = MatrixQ(0,1);
    cc[16] = MatrixQ(1,1);
    cc(16,0,0) = 2;
    cc[17] = MatrixQ(0,1);
    cc[18] = MatrixQ(1,1);
    cc(18,0,0) = 2;
  
    HomologyField ho = cc.homology();
    std::cout << ho << std::endl;
}

int main( int argc, char ** argv )
{
    std::cout.setf(std::ios::unitbuf);
    // test homology computation
    test_some_chain_complex_q();
    return 0;
    
    uint32_t errors = 0;
    // test whether the diagonalizer computes the rank correctly
    if( argc == 4 )
    {
        errors = test_rank_q( atoi(argv[1]), atoi(argv[2]), atoi(argv[3]) );
    }
    else if( argc >= 5 )
    {
        errors = test_rank_q( atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]) );
    }
    else
    {
        std::cerr << "Error: #Parameter < 3" << std::endl;
        return 1;
    }
    std::cout << "Total number of errors: " << errors << std::endl;
    
    return 0;
}

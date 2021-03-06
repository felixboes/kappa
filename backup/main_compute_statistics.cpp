#include <iostream>
#include <functional>
#include <string>

#include <cstddef> // has to be included before gmpxx with gcc 4.9
#include <gmpxx.h>

#include <libhomology/homology.hpp>

#include "kappa.hpp"

void print_usage(int, char** argv)
{
    std::cout << "Usage: " << argv[0] << " -g arg -m arg --first_diff arg --last_diff arg" << std::endl;
}

template< class MatrixT >
void print_statistics( MatrixT& M )
{

    // Watch out: Since we store the transpose of the differential,
    // we act as though rows were columns and vise versa.
    uint32_t num_rows = M.size1();
    uint32_t num_cols = M.size2();
    int32_t pos_infty = std::numeric_limits<int32_t>::max();
    int32_t neg_infty = std::numeric_limits<int32_t>::min();
    
    std::cout << "Number of rows: " << num_rows << std::endl;
    std::cout << "Number of columns: " << num_cols << std::endl;
    
    std::vector<mpz_class> num_entries_per_col( num_rows, 0 );
    std::vector<mpz_class> largest_entry_per_col( num_rows, neg_infty );
    std::vector<mpz_class> smallest_entry_per_col( num_rows, pos_infty );
    
    for( uint32_t i = 0; i < num_rows; ++i )
    {
        for( uint32_t j = 0; j < num_cols; ++j )
        {
            if( M(i,j) != 0 )
            {
                num_entries_per_col[i]++;
            }
            largest_entry_per_col[i] = std::max( largest_entry_per_col[i], mpz_class(M(i,j)) );
            smallest_entry_per_col[i] = std::min( smallest_entry_per_col[i], mpz_class(M(i,j)) );
        }
    }
    
    std::cout << "Number of non-zero entries (per column): |" << std::flush;
    for( auto& it : num_entries_per_col )
    {
        std::cout << it << "|";
    }
    std::cout << std::endl;
    
    std::cout << "Largest entries (per column): |";
    for( auto& it : largest_entry_per_col )
    {
        std::cout << it << "|";   
    }
    std::cout << std::endl;
    
    std::cout << "Smallest entries (per column): |";
    for( auto& it : smallest_entry_per_col )
    {
        std::cout << it << "|";
    }

    mpz_class total(0);
    for( uint32_t i = 0; i < num_rows; ++i )
    {
        total += num_entries_per_col[i];
    }
    if( num_rows > 0 )
    {
        std::cout << "Average number of non-zero entries (arithmetic mean): exact: " << mpq_class(total) / num_rows << " floored:" << total / num_rows << std::endl;
    }
    else
    {
        std::cout << "Average number of non-zero entries (arithmetic mean): exact: 0 floored: 0" << std::endl;
    }

}

int main(int argc, char** argv)
{
    // Parse configuration from command line arguments.
    SessionConfig conf(argc, argv);
    if( conf.option_set("help") )
    {
        print_usage(argc, argv);
        std::cout << conf.desc << std::endl;
        return 0;
    }
    
    if( ! ( conf.option_set( "gen" ) && conf.option_set( "pun" ) && conf.option_set( "first_diff" ) && conf.option_set( "last_diff" ) ) )
    {
        print_usage(argc, argv);
        std::cout << conf.desc << std::endl;
        return 1;
    }
    
    if ( conf.setup_configuration() == false )
    {
        std::cout << "The configuration could not been setup." << std::endl;
        return 2;
    }
    
    for( auto i =std::max<uint32_t>( conf.start_p, conf.num_punctures+1 ); i <= std::min<uint32_t>( conf.end_p, 4*conf.genus + 2*conf.num_punctures ); ++i )
    {
        #ifndef WE_USE_AN_OLD_COMPILER_THAT_DOES_NOT_SUPPORT_ALL_CPP_ELEVEN_FEATURES_OR_OPTIMIZATION
        MatrixZDontDiagonalize M = load_from_file_bz2<MatrixZDontDiagonalize>( "./cache/differentials/" + std::to_string(conf.genus) + "_" + std::to_string(conf.num_punctures) + "_" + std::to_string(i) );
        #else
        MatrixZDontDiagonalize M;
        load_from_file_bz2<MatrixZDontDiagonalize>( M, "./cache/differentials/" + std::string(conf.genus) + "_" + std::string(conf.num_punctures) + "_" + std::string(i) );
        #endif
        
        std::cout << "num rows " << M.size1() << " num cols " << M.size2() << std::endl;
        print_statistics(M);
        std::cout << std::endl;
//        std::cout << M << std::endl;
//        std::cout << std::endl;
        for (size_t i = 0; i < M.size1(); ++i)
        {
            for (size_t j = 0; j < M.size2(); ++j)
            {
                std::cout << (M(i,j) == 0 ? "0" : "*") << " ";
            }
            std::cout << std::endl;
        }
//        std::cout << M << std::endl;
//        BlockFinder<MatrixZDontDiagonalize> block_finder(M);
//        std::cout << "num blocks " << block_finder.num_blocks() << std::endl;
//        std::cout << "num non zero blocks " << block_finder.num_non_zero_blocks() << std::endl;
    }
    
    return 0;
}

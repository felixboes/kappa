#include <iostream>
#include <functional>
#include <string>

#include <gmpxx.h>

#include <homology.hpp>

#include "kappa.hpp"

void print_usage(int argc, char** argv)
{
    std::cout << "Usage: " << argv[0] << " genus num_punctures num_differential" << std::endl;
}

template< class MatrixT >
void print_statistics( MatrixT& M )
{
    uint32_t num_rows = M.size1();
    uint32_t num_cols = M.size2();
    int32_t pos_infty = std::numeric_limits<int32_t>::max();
    int32_t neg_infty = std::numeric_limits<int32_t>::min();
    
    std::cout << "Number of rows: " << num_rows << std::endl;
    std::cout << "Number of columns: " << num_cols << std::endl;
    
    std::vector<mpz_class> num_entries_per_col( num_cols, 0 );
    std::vector<mpz_class> largest_entry_per_col( num_cols, neg_infty );
    std::vector<mpz_class> smallest_entry_per_col( num_cols, pos_infty );
    
    for( uint32_t i = 0; i < num_cols; ++i )
    {
        for( uint32_t j = 0; j < num_rows; ++j )
        {
            if( M(j,i) != 0 )
            {
                num_entries_per_col[i]++;
            }
            largest_entry_per_col[i] = std::max( largest_entry_per_col[i], mpz_class(M(j,i)) );
            smallest_entry_per_col[i] = std::min( smallest_entry_per_col[i], mpz_class(M(j,i)) );
        }
    }
    
    std::cout << "Number of non-zero entries (per column): |";
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
    std::cout << std::endl;
    
    mpz_class total(0);
    for( uint32_t i = 0; i < num_cols; ++i )
    {
        total += num_entries_per_col[i];
    }
    if( num_cols > 0 )
    {
        std::cout << "Average number of non-zero entries (arithmetic mean): exact: " << mpq_class(total) / num_cols << " floored:" << total / num_cols << std::endl;
    }
    else
    {
        std::cout << "Average number of non-zero entries (arithmetic mean): exact: 0 floored: 0" << std::endl;
    }

}

int main(int argc, char** argv)
{
    if(argc < 4)
    {
        print_usage(argc, argv);
        return 1;
    }
    
    MatrixZDontDiagonalize M;
    load_from_file_bz2<MatrixZDontDiagonalize>( M, "./cache/differentials/" + std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) );
    std::cout << "num rows " << M.size1() << " num cols " << M.size2() << std::endl;
    /*print_statistics(M);

    for (int i = 0; i < M.size1(); ++i)
    {
        for (int j = 0; j < M.size2(); ++j)
        {
            std::cout << (M(i,j) == 0 ? "0" : "*") << " ";
        }
        std::cout << std::endl;
    }
    std::cout << M << std::endl;
    BlockFinder<MatrixZDontDiagonalize> block_finder(M);
    std::cout << "num blocks " << block_finder.num_blocks() << std::endl;
    std::cout << "num non zero blocks " << block_finder.num_non_zero_blocks() << std::endl;*/
    return 0;
}

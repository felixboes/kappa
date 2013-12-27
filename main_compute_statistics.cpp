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
    std::cout << "Number of rows: " << M.size1() << std::endl;
    std::cout << "Number of columns: " << M.size2() << std::endl;
    
    int32_t pos_infty = std::numeric_limits<int32_t>::max();
    int32_t neg_infty = std::numeric_limits<int32_t>::min();
    
    std::vector<mpz_class> num_entries_per_col(M.size2(), 0);
    std::vector<mpz_class>  largest_entry_per_col(M.size2(), neg_infty );
    std::vector<mpz_class> smallest_entry_per_col(M.size2(), pos_infty );
    
    for( uint32_t i = 0; i < M.size2(); ++i )
    {
        for( uint32_t j = 0; j < M.size1(); ++j )
        {
            if( M(j,i) != 0 )
            {
                num_entries_per_col[i]++;
                largest_entry_per_col[i] = std::max( largest_entry_per_col[i], mpz_class(M(j,i)) );
                smallest_entry_per_col[i] = std::min( smallest_entry_per_col[i], mpz_class(M(j,i)) );
            }
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
        if( it == neg_infty )
        {
            std::cout << "-inf|";
        }
        else
        {
            std::cout << it << "|";   
        }
    }
    std::cout << std::endl;
    
    std::cout << "Smallest entries (per column): |";
    for( auto& it : smallest_entry_per_col )
    {
        if( it == pos_infty )
        {
            std::cout << "+inf|";
        }
        else
        {
            std::cout << it << "|";   
        }
    }
    std::cout << std::endl;
    
    mpz_class total(0);
    for( uint32_t i = 0; i < M.size2(); ++i )
    {
        total += num_entries_per_col[i];
    }
    std::cout << "Average number of non-zero entries (arethmetic mean): exact: " << mpq_class(total) / M.size2() << " floored:" << total / M.size2() << std::endl;
}

int main(int argc, char** argv)
{
    if(argc < 4)
    {
        print_usage(argc, argv);
        return 1;
    }
    
    MatrixZDontDiagonalize M = load_from_file_bz2<MatrixZDontDiagonalize>( "./cache/differentials/" + std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) );
    std::cout << M << std::endl;
    print_statistics(M);
    
    return 0;
}

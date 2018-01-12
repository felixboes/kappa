#ifndef BLOCKFINDER_HPP
#define BLOCKFINDER_HPP

#include <stack>
#include <vector>

#include "tuple.hpp"

typedef std::vector< int > Component; ///< used to store a a subset of rows, or of columns of a matrix 
typedef std::pair< Component, Component > Block; ///< used to store a block of a matrix
typedef std::vector< Block > BlockPartition; ///< used to store a partition of the matrix into blocks

/**
 * @brief Determines the block structure of a given matrix of type MatrixT.
 * Thereby a block of a matrix is defined as a maximal collection C of rows and columns
 * of the matrix such that for each row r and each column c, there is an alternating 
 * sequence of rows and columns in C starting at r and ending at c such that for all consecutive pairs of 
 * a row r and a column c (or the other way around) the entry matrix(r, c) is non-zero. .
 */
template< class MatrixT >
class BlockFinder
{
public:
    typedef MatrixT MatrixType;
    
    /**
     * Determines a block partition of the given matrix, i.e. a collection of blocks
     * covering all rows of the matrix, and stores it in _block_part. 
     * @warning 0-rows build own blocks in the block partition.
     * @warning 0-columns are ignored in the block partition.
     */
    BlockFinder(MatrixType &);

    /**
     * @return _block_part
     */
    BlockPartition & block_partition();
    
    /**
     * @return Total number of blocks in _block_part.
     * @warning 0-rows build own blocks in the block partition.
     * @warning 0-columns are ignored in the block partition.
     */
    int num_blocks() const;

    /**
     * @return Number of blocks in _block_part which are not zero rows or zero blocks.
     */
    int num_non_zero_blocks() const;

protected:
    BlockPartition _block_part; ///< stores a block partition of the given matrix
};

#endif // BLOCKFINDER_HPP

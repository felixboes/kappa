#ifndef BLOCKFINDER_HPP
#define BLOCKFINDER_HPP

#include <stack>
#include <vector>

#include "tuple.hpp"

typedef std::vector< int > Component;
typedef std::pair< Component, Component > Block;
typedef std::vector< Block > BlockPartition;

template< class MatrixT >
class BlockFinder
{
public:
    typedef MatrixT MatrixType;
    /**
     * Determines a block partition of the given matrix and stores it in _block_part.
     * @warning 0-rows build own blocks in the block partition.
     */
    BlockFinder(MatrixType &);

    /**
     * @return _block_part
     */
    BlockPartition & block_partition();
    /**
     * @return Total number of blocks in _block_part.
     */
    int num_blocks();

    /**
     * @return Number of blocks in _block_part which are not zero rows.
     */
    int num_non_zero_blocks();

private:
    BlockPartition _block_part;
};

#include "blockfinder.ipp"

#endif // BLOCKFINDER_HPP

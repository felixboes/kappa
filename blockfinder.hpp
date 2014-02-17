#ifndef BLOCKFINDER_HPP
#define BLOCKFINDER_HPP

#include <stack>
#include <vector>

#include "tuple.hpp"

typedef std::vector< std::vector< Tuple > > BasisPartition;
typedef std::pair< BasisPartition, BasisPartition > Block;

template< class MatrixT >
class BlockFinder
{
public:
    typedef MatrixT MatrixType;

    BlockFinder();
    Block operator() ( MatrixType & );
    
};

#include "blockfinder.ipp"

#endif // BLOCKFINDER_HPP

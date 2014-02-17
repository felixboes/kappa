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

    BlockFinder();
    BlockPartition operator() ( MatrixType & );
};

#include "blockfinder.ipp"

#endif // BLOCKFINDER_HPP

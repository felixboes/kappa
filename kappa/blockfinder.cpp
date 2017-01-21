#include "blockfinder.hpp"
#include "blockfinder_impl.ipp"

#include "kappa.hpp"

template class BlockFinder< MatrixQ >;
template class BlockFinder< MatrixZm >;
template class BlockFinder< MatrixZDontDiagonalize >;

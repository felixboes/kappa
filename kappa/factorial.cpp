// The software kappa is a collection of programs to compute the homology of
// the moduli space of surfaces using the radial model.
// Copyright (C) 2013 - 2018  Felix Boes and Anna Hermann
// 
// This file is part of kappa.
// 
// kappa is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// kappa is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with kappa.  If not, see <http://www.gnu.org/licenses/>.


#include "factorial.hpp"

uint64_t factorial( const uint32_t k )
{
    switch( k )
    {
        case 0:  return 1LL;
        case 1:  return 1LL;
        case 2:  return 2LL;
        case 3:  return 6LL;
        case 4:  return 24LL;
        case 5:  return 120LL;
        case 6:  return 720LL;
        case 7:  return 5040LL;
        case 8:  return 40320LL;
        case 9:  return 362880LL;
        case 10: return 3628800LL;
        case 11: return 39916800LL;
        case 12: return 479001600LL;
        case 13: return 6227020800LL;
        case 14: return 87178291200LL;
        case 15: return 1307674368000LL;
        case 16: return 20922789888000LL;
        case 17: return 355687428096000LL;
        case 18: return 6402373705728000LL;
        case 19: return 121645100408832000LL;
        case 20: return 2432902008176640000LL;
        default : return 0; // 21! > 2^64 - 1
    }
}

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


#include "operationtester.hpp"
#include "operationtester_impl.ipp"

#define force_template_instantiation( CoefficientType, MatrixComplex, VectorType, TupleType  )\
    template class OperationTester< MatrixComplex,  VectorType  >;\
    template VectorType kappa_dual( const CoefficientType& c, const TupleType& t, const EhrBasis<TupleType>& b );\
    template void compute_and_add_kappa_dual_rec( const CoefficientType& c, const TupleType& t,const EhrBasis<TupleType>& b, VectorType& v, const std::vector<size_t> s, const size_t i );


force_template_instantiation(Q, ChainComplexQ, VectorQ, SymGrpTuple)
force_template_instantiation(Zm, ChainComplexZm, VectorZm, SymGrpTuple)

#undef force_template_instantiation

typedef OperationTesterCSS< ChainComplexQCSS,  VectorQ  > OperationTesterQCSS;
typedef OperationTesterCSS< ChainComplexZmCSS, VectorZm > OperationTesterZmCSS;

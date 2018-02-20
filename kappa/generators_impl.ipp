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


#include "generators.hpp"

template< class CoefficientT >
MonoCochainField< CoefficientT > create_cochain( const Generator& name )
{
    typedef MonoCochainField< CoefficientT > CochainType;
    switch( name )
    {
    
    case a:
    {
        CochainType cochain(0, 1, 2);
        cochain.set_name("a");
        cochain( create_cell(1, 2, 1) ) = CoefficientT(1);
        return cochain;
    }
    case b:
    {
        CochainType cochain(0, 2, 3);
        cochain.set_name("b");
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(2, 2, 1, 3, 1) );
        return cochain;
    }
    case c:
    {
        CochainType cochain(1, 0, 4);
        cochain.set_name("c");
        cochain( create_cell(2, 4, 2, 3, 1 ) ) = CoefficientT(1);
        return cochain;
    }
    case d:
    {
        CochainType cochain(1, 0, 3);
        cochain.set_name("d");
        cochain( create_cell(2, 3, 1, 2, 1 ) ) = CoefficientT(1);
        return cochain;
    }
    case e:
    {
        CochainType cochain(1, 1, 4);
        cochain.set_name("e");
        cochain.add_kappa_dual( CoefficientT(1),  create_cell(3, 3, 1, 4, 3, 2, 1) );
        return cochain;
    }
    case f:
    {
        CochainType cochain(1, 2, 5);
        cochain.set_name("f");
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 5, 4, 3, 1, 4, 3, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 5, 3, 4, 1, 5, 4, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 4, 1, 5, 4, 2, 1, 4, 3 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 3, 1, 5, 3, 2, 1, 5, 4 ) );
        
        return cochain;
    }
    case E_1_f:
    {
        CochainType cochain(2, 1, 6);
        cochain.set_name("E_1(f)");

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 5, 4, 6, 5, 3, 1, 5, 3, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 5, 3, 6, 5, 4, 1, 5, 3, 2, 1 ) );

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 5, 3, 6, 3, 4, 1, 6, 4, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 4, 3, 6, 3, 5, 1, 6, 4, 2, 1 ) );

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 5, 3, 4, 1, 6, 4, 2, 1, 4, 3 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 4, 3, 5, 1, 6, 4, 2, 1, 4, 3 ) );

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 5, 4, 3, 1, 6, 3, 2, 1, 5, 4 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 4, 3, 5, 1, 6, 3, 2, 1, 5, 4 ) );

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 5, 4, 3, 1, 6, 3, 2, 1, 6, 5 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 5, 3, 4, 1, 6, 3, 2, 1, 6, 5 ) );

        return cochain;
    }
    case E_2_f:
    {
        CochainType cochain(2, 1, 6);
        cochain.set_name("E_2(f)");

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 5, 3, 6, 4, 3, 1, 4, 3, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 5, 1, 6, 4, 3, 1, 4, 3, 2, 1 ) );

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 5, 4, 6, 3, 4, 1, 5, 4, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 5, 1, 6, 3, 4, 1, 5, 4, 2, 1 ) );

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 5, 4, 6, 3, 5, 1, 6, 5, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 4, 1, 6, 3, 5, 1, 6, 5, 2, 1 ) );

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 5, 4, 5, 1, 6, 5, 2, 1, 5, 3 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 4, 1, 5, 1, 6, 5, 2, 1, 5, 3 ) );

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 5, 3, 3, 1, 6, 3, 2, 1, 6, 4 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 5, 1, 3, 1, 6, 3, 2, 1, 6, 4 ) );

        return cochain;
    }
    case Qb:
    {
        CochainType cochain(0, 4, 5);
        cochain.set_name("Q(b)");
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 4, 3, 5, 3, 2, 1, 3, 1) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 4, 1, 5, 1, 3, 2, 4, 2) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 2, 1, 5, 1, 4, 3, 5, 3) );
        
        return cochain;
    }
    case Qc:
    {
        CochainType cochain(2, 0, 7);
        cochain.set_name("Q(c)");
        // Observe that every two adjecent cells differ by a face index 1, hence the sign is kept.
        cochain.add_kappa_dual( CoefficientT(1),  create_cell(4, 7, 5, 6, 4, 4, 2, 3, 1) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell(4, 7, 5, 6, 1, 4, 2, 3, 1) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell(4, 7, 5, 6, 2, 4, 2, 3, 1) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell(4, 7, 5, 6, 3, 4, 2, 3, 1) );
        // Now jumping to the groud differs by a face index of 4, hence we have to change the sign.
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 5, 6, 1, 5, 3, 4, 2) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 2, 6, 1, 5, 3, 4, 2) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 3, 6, 1, 5, 3, 4, 2) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 4, 6, 1, 5, 3, 4, 2) );
        // Now jumping to the groud differs by a face index of 4, hence we have to change the sign.
        cochain.add_kappa_dual( CoefficientT(1),  create_cell(4, 7, 2, 6, 1, 6, 4, 5, 3) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell(4, 7, 2, 3, 1, 6, 4, 5, 3) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell(4, 7, 2, 4, 1, 6, 4, 5, 3) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell(4, 7, 2, 5, 1, 6, 4, 5, 3) );
        // Now jumping to the groud differs by a face index of 4, hence we have to change the sign.
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 2, 3, 1, 7, 5, 6, 4) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 4, 2, 3, 1, 7, 5, 6, 4) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 5, 2, 3, 1, 7, 5, 6, 4) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 6, 2, 3, 1, 7, 5, 6, 4) );
        
        return cochain;
    }
    case Qd:
    {
        CochainType cochain(2, 0, 5);
        cochain.set_name("Q(d)");
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 3, 3, 1, 2, 1 ) );
    
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 1, 3, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 4, 1, 3, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 2, 3, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 4, 2, 3, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 2, 4, 2, 3, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 4, 1, 4, 2, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 2, 4, 1, 4, 2, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 1, 4, 2, 3, 2 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 2, 1, 4, 2, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 2, 1, 4, 2, 3, 2 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 3, 1, 4, 2, 3, 2 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 2, 1, 5, 3, 4, 3 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 3, 1, 2, 1, 5, 3, 4, 3 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 4, 1, 2, 1, 5, 3, 4, 3 ) );
        
        return cochain;
    }
    case Te:
    {
        CochainType cochain(2, 0, 5);
        cochain.set_name("T(e)");
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 5, 3, 3, 1, 4, 3, 2, 1) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 5, 1, 3, 1, 4, 3, 2, 1) );
        return cochain;
    }
    case Eb:
    {
        CochainType cochain(1, 1, 3);
        cochain.set_name("E(b)");
        cochain.add_kappa_dual( CoefficientT(1), create_cell(3, 2, 1, 2, 1, 3, 1) );
        return cochain;
    }
    case TEb:
    {
        CochainType cochain(2, 0, 4);
        cochain.set_name("T(E(b))");
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 4, 2, 2, 1, 2, 1, 3, 1) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 4, 1, 2, 1, 2, 1, 3, 1) );
        return cochain;
    }
    case Q_alpha_inv_c:
    {
        CochainType cochain( create_cochain<CoefficientT>( Generator::Qc ) );
        cochain.set_name("Q_alpha^{-1}(c)");
        return cochain;
    }
    case Q_beta_c:
    {
        CochainType cochain( 2, 0, 7 );
        cochain.set_name("Q_beta(c)");
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell(4, 5, 2, 3, 1, 7, 5, 6, 4) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell(4, 4, 2, 3, 1, 7, 5, 6, 4) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell(4, 7, 2, 3, 1, 7, 5, 6, 4) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 2, 5, 1, 6, 4, 5, 3) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 4, 5, 1, 6, 3, 4, 2) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 4, 5, 3, 6, 2, 3, 1) );
        
        return cochain;
    }
    case Q_gamma_c:
    {
        CochainType cochain( 2, 0, 7 );
        cochain.set_name("Q_gamma(c)");
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell(4, 5, 2, 6, 4, 7, 2, 3, 1) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell(4, 5, 1, 6, 4, 7, 2, 3, 1) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell(4, 7, 5, 6, 4, 7, 2, 3, 1) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 4, 5, 3, 6, 2, 3, 1) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 4, 5, 1, 6, 3, 4, 2) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell(4, 7, 2, 5, 1, 6, 4, 5, 3) );
        
        return cochain;
    }
    case Q_alpha_inv_d:
    {
        CochainType cochain( create_cochain<CoefficientT>( Generator::Qd ) );
        cochain.set_name("Q_alpha^{-1}(d)");
        return cochain;
    }
    case Q_beta_d:
    {
        CochainType cochain(2, 0, 5);
        cochain.set_name("Q_beta(d)");
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 3, 1, 2, 1, 5, 3, 4, 3 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 2, 1, 5, 3, 4, 3 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 3, 1, 4, 2, 3, 2 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 2, 3, 2, 4, 1, 2, 1 ) );
        
        return cochain;
    }
    case Q_gamma_d:
    {
        CochainType cochain(2, 0, 5);
        cochain.set_name("Q_gamma(d)");
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 3, 1, 4, 3, 5, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 3, 5, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 2, 3, 2, 4, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 3, 1, 4, 2, 3, 2 ) );
        
        return cochain;
    }
    case R_a_e:
    {
        CochainType cochain(1, 2, 5);
        cochain.set_name("R(a,e)");
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 5, 4, 3, 1, 4, 3, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 5, 2, 3, 1, 4, 3, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 5, 1, 4, 2, 5, 4, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 3, 1, 4, 2, 5, 4, 3, 2 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 4, 2, 5, 4, 3, 2, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 4, 1, 5, 4, 3, 1, 3, 2 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 4, 1, 5, 4, 2, 1, 4, 3 ) );
        
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 3, 1, 5, 3, 2, 1, 5, 4 ) );
        
        return cochain;
    }
    case R_alpha_inv_beta_c_d:
    {
        CochainType cochain(2, 0, 6);
        cochain.set_name("R(alpha^{-1}*beta, c, d)");
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 5, 3, 3, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 1, 3, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 3, 4, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 1, 5, 3, 4, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 3, 5, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 3, 1, 5, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 5, 1, 2, 1, 6, 4, 5, 3 ) );
        
        return cochain;
    }
    case R_alpha_inv_gamma_inv_c_d:
    {
        CochainType cochain(2, 0, 6);
        cochain.set_name("R(alpha^{-1}*gamma^{-1}, c, d)");
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 5, 3, 3, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 1, 3, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 3, 4, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 1, 5, 3, 4, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 3, 5, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 3, 1, 5, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 5, 3, 6, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 4, 1, 5, 3, 6, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 1, 2, 1, 6, 4, 5, 3 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 3, 1, 2, 1, 6, 4, 5, 3 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 4, 1, 2, 1, 6, 4, 5, 3 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 5, 1, 2, 1, 6, 4, 5, 3 ) );
        
        return cochain;
    }
    case R_alpha_inv_gamma_c_d:
    {
        CochainType cochain(2, 0, 6);
        cochain.set_name("R(alpha^{-1}*gamma, c, d)");
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 5, 3, 3, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 1, 3, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 3, 4, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 1, 5, 3, 4, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 3, 5, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 3, 1, 5, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 3, 4, 2, 5, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 3, 4, 1, 5, 2, 3, 2 ) );
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 2, 4, 1, 5, 3, 4, 3 ) );
        
        // horizontal swap
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 5, 3, 4, 3, 6, 2, 3, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 5, 1, 4, 1, 6, 3, 4, 2 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 5, 1, 2, 1, 6, 4, 5, 3 ) );        
        
        return cochain;
    }
    case R_alpha_inv_beta_inv_c_d:
    {
        CochainType cochain(2, 0, 6);
        cochain.set_name("R(alpha^{-1}*beta^{-1}, c, d)");
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 5, 3, 3, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 1, 3, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 3, 4, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 1, 5, 3, 4, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 3, 5, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 3, 1, 5, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 5, 3, 6, 1, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 4, 1, 5, 3, 6, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 5, 3, 4, 2, 6, 1, 2, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 5, 3, 4, 1, 6, 2, 3, 2 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 5, 2, 4, 1, 6, 3, 4, 3 ) );
        
        // horizontal swap
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 3, 4, 3, 5, 2, 3, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 1, 4, 1, 5, 3, 4, 2 ) );
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 1, 2, 1, 5, 4, 5, 3 ) );  
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 1, 2, 1, 6, 4, 5, 3 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 3, 1, 2, 1, 6, 4, 5, 3 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 4, 1, 2, 1, 6, 4, 5, 3 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 5, 1, 2, 1, 6, 4, 5, 3 ) );
        
        return cochain;
    }
    case R_alpha_inv_alpha_inv_c_d:
    {
        CochainType cochain(2, 0, 6);
        cochain.set_name("R(alpha^{-1}*alpha^{-1}, c, d)");
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 5, 3, 3, 1, 2, 1 ) ); // 3
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 1, 3, 1, 2, 1 ) ); // 1
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 3, 4, 1, 2, 1 ) ); // 4
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 1, 5, 3, 4, 1, 2, 1 ) ); // 1
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 3, 5, 1, 2, 1 ) ); // 5
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 3, 1, 5, 1, 2, 1 ) ); // 1
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 5, 3, 6, 1, 2, 1 ) ); // 6
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 4, 1, 5, 3, 6, 1, 2, 1 ) ); // 1
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 5, 3, 4, 2, 6, 1, 2, 1 ) ); // 3
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 5, 3, 4, 1, 6, 2, 3, 2 ) ); // 1 <-> 3
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 5, 2, 4, 1, 6, 3, 4, 3 ) ); // 2 <-> 4
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 5, 2, 3, 1, 6, 4, 5, 4 ) ); // 3 <-> 5
        
        // horizontal swap (number of sign changes is even)
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 4, 4, 2, 3, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 1, 4, 2, 3, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 1, 5, 1, 4, 2, 3, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 2, 4, 2, 3, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 1, 5, 2, 4, 2, 3, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 2, 5, 2, 4, 2, 3, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 4, 5, 3, 4, 2, 3, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 1, 5, 3, 4, 2, 3, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 2, 5, 3, 4, 2, 3, 1 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 3, 5, 3, 4, 2, 3, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 2, 5, 1, 5, 3, 4, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 3, 5, 1, 5, 3, 4, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 5, 1, 5, 3, 4, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 1, 5, 1, 5, 3, 4, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 3, 2, 1, 5, 3, 4, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 2, 1, 5, 3, 4, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 1, 2, 1, 5, 3, 4, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 4, 3, 1, 5, 3, 4, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 1, 3, 1, 5, 3, 4, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 6, 1, 4, 1, 5, 3, 4, 2 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 1, 2, 1, 6, 4, 5, 3 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 3, 1, 2, 1, 6, 4, 5, 3 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 4, 1, 2, 1, 6, 4, 5, 3 ) );
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 5, 1, 2, 1, 6, 4, 5, 3 ) );
        
        return cochain;
    }
    case R_alpha_inv_beta_c_c:
    {
        CochainType cochain(2, 0, 7);
        cochain.set_name("R(alpha^{-1}*beta, c, c)");
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 7, 5, 6, 4, 4, 2, 3, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 7, 5, 6, 1, 4, 2, 3, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 7, 5, 6, 2, 4, 2, 3, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 7, 5, 6, 4, 5, 2, 3, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 7, 1, 6, 4, 5, 2, 3, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 7, 2, 6, 4, 5, 2, 3, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 7, 5, 6, 4, 6, 2, 3, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 7, 5, 4, 1, 6, 2, 3, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1),  create_cell( 4, 7, 5, 4, 2, 6, 2, 3, 1 ) );
        
        cochain.add_kappa_dual( CoefficientT(-1), create_cell( 4, 6, 2, 3, 1, 7, 5, 6, 4 ) );
        
        return cochain;
    }
    case R_a_Te:
    {
        CochainType cochain(2, 1, 6);
        cochain.set_name("R(a,T(e))");

        // Fiddle first slit of 'a' down in first half of 'T(e)'
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 5, 5, 3, 3, 1, 4, 3, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 1, 5, 3, 3, 1, 4, 3, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 3, 5, 3, 3, 1, 4, 3, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 4, 5, 3, 3, 1, 4, 3, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 2, 5, 3, 3, 1, 4, 3, 2, 1 ) );

        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 1, 6, 4, 4, 2, 5, 4, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 2, 1, 6, 4, 4, 2, 5, 4, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 4, 1, 6, 4, 4, 2, 5, 4, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 5, 1, 6, 4, 4, 2, 5, 4, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 3, 1, 6, 4, 4, 2, 5, 4, 3, 2 ) );

        // Fiddle second slit of 'a' down in first half of 'T(e)'
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 5, 5, 1, 3, 1, 4, 3, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 3, 5, 1, 3, 1, 4, 3, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 1, 5, 1, 3, 1, 4, 3, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 4, 5, 1, 3, 1, 4, 3, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 2, 5, 1, 3, 1, 4, 3, 2, 1 ) );

        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 1, 6, 2, 4, 2, 5, 4, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 4, 1, 6, 2, 4, 2, 5, 4, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 2, 1, 6, 2, 4, 2, 5, 4, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 5, 1, 6, 2, 4, 2, 5, 4, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 3, 1, 6, 2, 4, 2, 5, 4, 3, 2 ) );

        // Now vertical movement.

        // Fiddle first half of 'T(e) through a.
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 4, 4, 2, 5, 4, 3, 2, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 4, 4, 1, 5, 4, 3, 1, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 4, 4, 1, 5, 4, 2, 1, 4, 3 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 3, 3, 1, 5, 3, 2, 1, 5, 4 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 3, 3, 1, 4, 3, 2, 1, 6, 5 ) );

        // Fiddle second half of 'T(e) through a.
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 2, 4, 2, 5, 4, 3, 2, 2, 1 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 1, 4, 1, 5, 4, 3, 1, 3, 2 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 1, 4, 1, 5, 4, 2, 1, 4, 3 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 1, 3, 1, 5, 3, 2, 1, 5, 4 ) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell( 5, 6, 1, 3, 1, 4, 3, 2, 1, 6, 5 ) );

        // Now vertical movement.

        return cochain;
    }
    case T_1_f:
    {
        CochainType cochain(2, 1, 6, false, "T_1(f)");

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 6, 4, 5, 4, 3, 1, 4, 3, 2, 1) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 6, 3, 5, 3, 4, 1, 5, 4, 2, 1) );

        // Now vertical movement.

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 6, 3, 4, 1, 5, 4, 2, 1, 4, 3) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 6, 4, 3, 1, 5, 3, 2, 1, 5, 4) );

        return cochain;
    }
    case T_2_f:
    {
        CochainType cochain(2, 1, 6, false, "T_2(f)");

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 6, 3, 5, 4, 3, 1, 4, 3, 2, 1) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 6, 4, 5, 3, 4, 1, 5, 4, 2, 1) );

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 6, 1, 5, 4, 3, 1, 4, 3, 2, 1) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 6, 1, 5, 3, 4, 1, 5, 4, 2, 1) );

        // Now vertical movement.

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 6, 4, 4, 1, 5, 4, 2, 1, 4, 3) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 6, 3, 3, 1, 5, 3, 2, 1, 5, 4) );

        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 6, 1, 4, 1, 5, 4, 2, 1, 4, 3) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(5, 6, 1, 3, 1, 5, 3, 2, 1, 5, 4) );

        return cochain;
    }
    case radial_a:
    {
        CochainType cochain(0, 2, 2, true, "a");

        cochain.add_kappa_dual( CoefficientT(1), create_cell(1, 2, 1) );

        return cochain;
    }
    case radial_b:
    {
        CochainType cochain(0, 2, 1, true, "b");

        cochain.add_kappa_dual( CoefficientT(1), create_cell(1, 1, 0) );

        return cochain;
    }
    case radial_Srad_Te:
    {
        CochainType cochain(2, 1, 4, true, "S(rad(T(e))");

        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 4, 0, 2, 0, 3, 2, 1, 0) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 4, 3, 4, 1, 2, 1, 4, 0) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 3, 2, 3, 0, 1, 0, 4, 3) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 2, 1, 4, 2, 4, 0, 3, 2) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 1, 0, 3, 1, 4, 3, 2, 1) );

        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 4, 2, 2, 0, 3, 2, 1, 0) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 3, 1, 4, 1, 2, 1, 4, 0) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 2, 0, 3, 0, 1, 0, 4, 3) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 4, 1, 4, 2, 4, 0, 3, 2) );
        cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 3, 0, 3, 1, 4, 3, 2, 1) );

        return cochain;
    }
    
    }
    
    
    return CochainType(0,1,2, false, "No Cochain loaded. Did you forget the return statement in generators_impl.hpp?");
}



# The software pykappa is a bunch of programs to compute the homology of
# Sullivan diagrams.
# Copyright (C) 2015, 2016  Felix Boes
#
# This file is part of pykappa.
#
# pyradbar is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyradbar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pykappa.  If not, see <http://www.gnu.org/licenses/>.

from perm import *


class Cell:

    def __init__(self, p, g, m):
        self._p = p
        self._g = g
        self._m = m
        self._h = 2*self._g + self._m
        self._inhomogenous = h*[ Permutation(p) ]


    def __getitem__(self, j):
        return self._inhomogenous[j+1]


    def norm(self):
        return sum( [ self._inhomogenous[i].norm() for i in range(self._h) ] )


    def num_cycles(self):
        sigma_h = Permutation(self._p)

        for tau in reversed( self._inhomogenous ):
            sigma_h *= tau
        sigma_h *= Permutation.get_long_cycle(self._p)

        return sigma_h.num_cyc()


    def has_correct_num_cycles(self):
        return self.num_cycles() == self._m


    def connected_components(self):
        return


    def num_clusters(self):
        return


    def monotone(self):
        return


    def f(self, i):
        return


    def phi(self, q, i):
        return


    def d_hor(self, i):
        return


    def d_hor_reduced(self, i):
        return


    def d_hor_double_complex(self, i):
        return


    def orientation_sign(self):
        return


    def __eq__(self, other):
        return


    def __bool__(self):
        return

    # In Python 2, the evalutation to bool is done by the function __nonzero__.
    __nonzero__ = __bool__


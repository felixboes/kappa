# The software pykappa is a bunch of programs to compute the homology of
# Sullivan diagrams.
# Copyright (C) 2015, 2016  Felix Boes
#
# This file is part of pykappa.
#
# pykappa is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pykappa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pykappa.  If not, see <http://www.gnu.org/licenses/>.

class Transposition:

    def __init__(self, a, b):
        self._a = max(a, b)
        self._b = min(a, b)

    def get_a(self):
        return self._a

    def get_b(self):
        return self._b

    def __str__(self):
        return '(' + str(self._a) + ',' + str(self._b) + ')'

class Permutation:

    # Static variable that is often used.
    long_cycle = {}
    long_cycle_inv = {}

    @classmethod
    def get_long_cycle(cls, p):
        if not cls.long_cycle.has_key(p):
            long_cyc_p = Permutation(p)
            long_cyc_p._repr = [i + 1 for i in range(p)] + [0]
            cls.long_cycle[p] = long_cyc_p
        return cls.long_cycle[p]

    @classmethod
    def get_long_cycle_inv(cls, p):
        if not cls.long_cycle_inv.has_key(p):
            long_cyc_inv_p = Permutation(p)
            long_cyc_inv_p._repr = [i + 1 for i in range(p)] + [0]
            cls.long_cycle_inv[p] = long_cyc_inv_p
        return cls.long_cycle_inv[p]

    def __init__(self, p, tr=None):
       self._p = p
       self._repr = [ i for i in range(self._p + 1) ]
       if tr is not None:
           self._repr[tr.get_a()] = tr.get_b()
           self._repr[tr.get_b()] = tr.get_a()

       self._valid = True

    def __mul__(self, other):
        if self._p != other._p:
            res = Permutation(0)
            res._valid = False
            return res

        res = Permutation(self._p)
        res._repr = [ self._repr[ other._repr[ i ] ] for i in range(self._p + 1) ]
        return res

    def norm(self):
        i = 0
        norm = self._p+1
        visited = (self._p + 2) * [False]
        visited[self._p + 1] = True

        while visited[i] == False:
            norm -= 1
            j = self._repr[i]
            visited[j] = True
            while j != i:
                j = self._repr[j]
                visited[j] = True

            while visited[i] == True and i <= self._p:
                i += 1
        return norm


    def num_cyc(self):
        return (self._p+1) - self.norm()


    def fixed_pts(self):
        return [x for x in range(self._p+1) if self._repr[x] == x]


    def num_fixed_pts(self):
        return len( self.fixed_pts() )


    def cycle_decomposition(self):
        i = 0
        cycle_decomp = []
        visited = (self._p + 2) * [False]
        visited[self._p + 1] = True

        while visited[i] == False:
            cycle_decomp.append([i])
            j = self._repr[i]
            visited[j] = True
            while j != i:
                cycle_decomp[-1].append(j)
                j = self._repr[j]
                visited[j] = True

            while visited[i] == True and i <= self._p:
                i += 1

        return tuple(tuple(x) for x in cycle_decomp)


    def d(self, num_i):
        if self._p > 1:
            i = num_i % (self._p+1)

            suc = self._repr[i]
            transp = [x for x in range(self._p+1)]
            transp[i] = suc
            transp[suc] = i
            intermediate_perm = [ transp[self._repr[x]] for x in range(self._p+1) if x != i ]

            self._repr = [x if x < i else x - 1 for x in intermediate_perm]
            self._p -= 1
        else:
            self._valid = False


    def __eq__(self, other):
        if self._valid == other._valid == True and self._p == other._p and self._repr == other._repr:
            return True
        else:
            return False


    def __str__(self):
        return str(self.cycle_decomposition())
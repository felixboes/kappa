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


import copy
import gzip
import inspect
import networkx as nx
import os
import cPickle as pickle
import sys

from perm import *

class RadialCell:
    def __init__(self, h, tau):
        self._p = 2*h
        self._h = h
        self._inhomogenous = [ Transposition(tau_i[0], tau_i[1]) for tau_i in tau ]
        self._m = None
        self._g = None
        self._hash = None

        self.compute_g_and_m()

        self._valid = True

    def __getitem__(self, j):
        return self._inhomogenous[j-1]

    def __setitem__(self, j, val):
        self._inhomogenous[j-1] = val

    def degree(self):
        return self._p

    def get_g(self):
        return self._g

    def get_m(self):
        return self._m

    def get_h(self):
        return self._h

    def norm(self):
        return sum( [ Permutation(self._p, self._inhomogenous[i]).norm() for i in range(self._h) ] )

    def compute_m(self):
        sigma_h = Permutation(self._p)

        for tau_i in reversed( self._inhomogenous ):
            t = Permutation(self._p, tau_i)
            sigma_h = sigma_h * t
        sigma_h = sigma_h * Permutation.get_long_cycle(self._p)

        return sigma_h.num_cyc()

    def compute_g_and_m(self):
        self._m = self.compute_m()
        self._g = (self._h - self._m + 1) / 2

    def has_correct_num_cycles(self):
        return self.compute_m() == self._m

    def connected_components(self):
        G = nx.Graph()
        for edge in self._inhomogenous:
            G.add_edge(edge.get_a(), edge.get_b())

        return nx.connected_components(G)

    def num_clusters(self):
        G = nx.Graph()
        for edge in self._inhomogenous:
            G.add_edge( edge.get_a(), edge.get_b() )

        return nx.number_connected_components(G)

    def monotone(self):
        for q in range( 1, self._h, 1 ):
            if self[q+1].get_a() < self[q].get_a():
                return False
        return True

    def f(self, i):
        # We denote tau_{i+1} | tau_i by (ab)(cd).

        a,b = self[i+1].get_a(), self[i+1].get_b()
        c,d = self[i].get_a(), self[i].get_b()

        if a == c: # (a )(a )
            if b == d: # (ab)(ab) = id and (a*)(a*) = 0
                self._inhomogenous = None
                self._valid = False
                return False
            elif a == b: # (a*)(ad) = (ad*) -> (d*)(ad)
                self[i+1] = Transposition(d, d)
                return True
            elif a == d: # (ab)(a*) = (ab*) -> (b*)(ab)
                self[i+1] = Transposition(b, b)
                self[i]   = Transposition(a, b)
                return True
            else: # (ab)(ad) = (adb) = (bd)(ab) or (db)(ab)
                self[i+1] = Transposition(b, d)
                self[i]   = Transposition(a, b)
                return True
        elif a == b: # (a*)(  )
            if a == d: # (a*)(ca) = (ca*) -> (a*)(ca)
                return True
            else: # (a*)(c*) or (a*)(cd) and the two transpositions commute
                if( a < c ):
                    return True
                else:
                    self[i], self[i+1] = self[i+1], self[i]
                    return True
        # Case: a > b
        elif b == c: # (ab)(b )
            if c == d: # (ab)(b*) = (ab*) -> (b*)(ab)
                self[i], self[i+1] = self[i+1], self[i]
                return True
            else: # (ab)(bd) = (abd) and (abd)' = (ad) and (abd)(ad) = (bd)
                self[i+1] = Transposition(b, d)
                self[i]   = Transposition(a, d)
                return True
            # These are all cases since a > b >= d.
        elif a == d: # (ab)(ca) = (acb) and c > a -> (ab)(ca)
            return True
        elif b == d: # (ab)(cb) = (abc) -> a < c: (ab)(cb) oder a > c:(cb)(ac)
            if a < c:
                return True
            else:    # (abc)' = (ac) and (abc)(ac) = (cb)
                self[i+1] = Transposition(c, b)
                self[i]   = Transposition(a, c)
                return True
        else: # (ab)(c ) with d not a and not b, hence disjoint transpositions
            if a < c:
                return True
            else:
                self[i], self[i+1] = self[i+1], self[i]
                return True

    def phi(self, q, i):
        if i == 0 or i > q:
            raise RuntimeError('q=' + str(q) + ' and i=' + str(i) + ', but 0 < i <= q is asserted.')

        for j in range( q-1, i-1, -1):
            if self.f(j) == False:
                return False

        return True

    def d_hor(self, k):
        sigma = Permutation.get_long_cycle(self._p)._repr
        sigma_inv = Permutation.get_long_cycle_inv(self._p)._repr

        for q in range(1, self._h+1):
            # We denote tau_q = (a,b)
            a, b = self[q].get_a(), self[q].get_b()

            # We denote (k, sigma_{q-1}(k)) = (k,l)
            l = sigma[k]

            # Compute tau':
            # Most of the time the transpositions are disjoint hence (a,b)(k,l) = (k,l)(a,b) and
            # the left transposition will be killed by D_k
            if k != a and k != b and l != a and l != b:
                pass
            # The degenerate case:
            # k and l are both part of tau_q.
            elif a == max(k, l) and ( b == min(k, l) or k == l):
                return False
            # The non degenerate case:
            else:
                # Compute Z = (a,b)(k,l)
                # Z(k) = l iff l != a and l != b hence k == a or k == b
                # In this case: (Z(k),k)Z = (k,l)(a,b)(k,l) = (c,l) with c != k
                if l != a and l != b:
                    if a != k:
                        self[q] = Transposition(a, l)
                    else:
                        self[q] = Transposition(b, l)
                # Z(k) != l iff l = a or l = b hence k != a and k != b
                # In this case: (Z(k),k)Z = (c,k)(a,b)(k,l) = (c,l) with c != l,
                # but we see that this is just (a,b):
                # if l != a, then k,l != a and we map
                # a to a to b != c since otherwise a would map to k
                # therefore a = c and l = b.
                # Thus we do not need to alter the boundary.

            # Compute sigma_{q}
            # (a,b)sigma(k) =
            #   a          if k = sigma^{-1}(b)
            #   b          if k = sigma^{-1}(a)
            #   sigma(k)   else
            # this is done by swapping the values of sigma^{-1}(a) and sigma^{-1}(b) under sigma:
            sigma[ sigma_inv[a] ], sigma[ sigma_inv[b] ] = sigma[ sigma_inv[b] ], sigma[ sigma_inv[a] ]

            # sigma^{-1}(a,b) (k) =
            #   sigma^{-1}(b)  if k = a
            #   sigma^{-1}(a)  if k = b
            #   sigma^{-1}(k)  else
            # this is done by interchanging the values of a and b under sigma^{-1}
            sigma_inv[a], sigma_inv[b] = sigma_inv[b], sigma_inv[a]


        #boundary is monotone iff its renormalization is monotone. Thus we check for monotony now to avoid unnecessary renormalization.
        if not self.monotone():
            return False

        # Renormalize all tau'
        for q in range(1, self._h+1):
            # We denote tau_q = (a,b)
            a, b = self[q].get_a(), self[q].get_b()

            if a > k:
                a -= 1
            if b > k:
                b -= 1
            self[q] = Transposition(a,b)

        self._p -= 1

        return True

    def face(self, i):
        raise NotImplementedError('Complete the code that is commented out.')

    def sigma_h(self):
        # initialize with sigma_0
        sigma_inv = Permutation.get_long_cycle_inv(self._p)

        for i in range(1, self._h+1):
            # We denote tau_i = (a, b)
            a, b = self[i].get_a(), self[i].get_b()

            # compute sigma_i^{-1}
            sigma_inv._repr[ a ], sigma_inv._repr[ b ] = sigma_inv._repr[ b ], sigma_inv._repr[ a ]

        # compute sigma
        sigma = Permutation(self._p)
        for i in range(0, self._p+1):
            sigma._repr[sigma_inv._repr[i]] = i

        return sigma

    def orientation_sign(self):
        sigma = self.sigma_h()
        cycle_decomp = sigma.cycle_decomposition()
        num_cycles = len(cycle_decomp)
        sign = {}

        i = 1 # counter of cycles
        for it, cycle in enumerate(cycle_decomp):
            min_symbol = cycle[0]
            # if the cycle is a fixpoint (a), we set sign(a) = 0 for the sake of completeness.
            if len(cycle) == 1:
                sign[min_symbol] = 0
                i += 1
                continue

            # determine the second min symbol of the cycle
            second_min_symbol = min( cycle[1:] )

            # Find k.
            # note that
            #   a_{1,1} < ... < a_{i,1} < b ,
            # hence
            #   a_{i-l,1} < b    for    l >= 0
            # is impossible and we can start to search at the position k = i.
            k = copy.deepcopy(i)

            # note that since we exclude the case that second_min_sybols will be sorted in at the end,
            # min_symbol can be sorted in between the cycles and this loop will never reach the liniel.

            # note that initially (for k = i) it_2.first = a_{k,1} < b and we found our position k iff the first time b < a_{k+1,1}
            # in the other case we have again a_{k+1,1} < b. Induction.
            for it_2 in range(it+1, num_cycles):
                next_min_symbol = cycle_decomp[it_2][0]
                if second_min_symbol < next_min_symbol:
                    if ((k - i) % 2) == 0:
                        sign[min_symbol] = 1
                    else:
                        sign[min_symbol] = -1
                    break
                k += 1

            # If b > a_{m, 1}, we still need to determine the sign of min_symbol.
            if k == num_cycles:
                if ((k - i) % 2) == 0:
                    sign[min_symbol] = 1
                else:
                    sign[min_symbol] = -1

            # for all other symbols of the cycle, we set the sign to 1
            for c in cycle:
                if c != min_symbol:
                    sign[c] = 1
            i += 1
        return sign

    def get_clean_copy(self):
        ret = copy.deepcopy(self)
        ret._hash = None
        return ret

    def __hash__(self):
        if self._hash is None:
            self._hash = int(0)
            offset = int(2)
            for tau_i in self._inhomogenous:
                self._hash += offset * ( tau_i.get_a() + 8 * tau_i.get_b() )
                offset *= 16
            return self._hash
        else:
            hash_tmp = int(0)
            offset = int(2)
            for tau_i in self._inhomogenous:
                hash_tmp += offset * ( tau_i.get_a() + 8 * tau_i.get_b() )
                offset *= 16
            if self._hash == hash_tmp:
                return self._hash
            else:
                raise ValueError("Hash values are different..." + str(self._hash) + " " + str(hash_tmp) + "\n")

    def __eq__(self, other):
        if self._valid != other._valid:
            return False
        else:
            return (
               self._p == other._p and
               self._g == other._g and
               self._m == other._m and
               self._inhomogenous == other._inhomogenous
            )

    def __bool__(self):
        return self._valid

    # In Python 2, the evalutation to bool is done by the function __nonzero__.
    __nonzero__ = __bool__

    def __str__(self):
        s = ''
        for q in range(self._h, 1, -1):
          s += str(self[q]) + '|'
        s += str(self[1])

        return s


class TopRadialCellTauArchive:
    def __init__(self, gziparchive=None):
        self._archive = gziparchive

    def __iter__(self):
        return self

    def next(self):
        if self._archive is None:
            raise StopIteration

        try:
            tau = pickle.load(self._archive)
        except EOFError:
            # We reached the archives end.
            raise StopIteration
        except:
            # Other error.
            frameinfo = inspect.getframeinfo(inspect.currentframe())
            sys.stdout.write(str(frameinfo.filename) + ' ' + str(frameinfo.lineno) + '\n')
            e, p, t = sys.exc_info()
            sys.stdout.write(str(e) + ' ' + str(p) + '\n')
            raise StopIteration
        else:
            return tau

    # Python 3:
    __next__ = next

    def close(self):
        if isinstance(self._archive, gzip.GzipFile):
            self._archive.close()
        elif self._archive is not None:
            sys.stdout.write("Wrong file type: " + str(type(self._archive)) + '\n')


class LoadTopRadialCellTau:
    def __init__(self, g, m):
        self._g = g
        self._m = m
        self._valid = True
        self._archive = None
        if g < 0 or m < 1:
            self._valid = False

    def __enter__(self):
        # open archive
        try:
            if (os.path.exists('./data/') and os.path.exists('./data/top_cell_g_' + str(self._g) + '_m_' + str(self._m) + '.bz2')):
                self._archive = TopRadialCellTauArchive(gzip.GzipFile('./data/top_cell_g_' + str(self._g) + '_m_' + str(self._m) + '.bz2','rb'))
                return self._archive
            else:
                sys.stdout.write('File ' + './data/top_cell_g_' + str(self._g) + '_m_' + str(self._m) + '.bz2 does not name a valid file.\n')
                self._valid = False
                return None
        except:
            # Print the error.
            frameinfo = inspect.getframeinfo(inspect.currentframe())
            sys.stdout.write(str(frameinfo.filename) + ' ' + str(frameinfo.lineno) + '\n')
            e, p, t = sys.exc_info()
            sys.stdout.write(str(e) + ' ' + str(p) + '\n')
            return None

    def __exit__(self, exc_type, exc_val, exc_tb):

        if exc_type is not None:
            if isinstance(self._archive, TopRadialCellTauArchive):
                self._archive.close()
            return False
        else:
            if isinstance(self._archive, TopRadialCellTauArchive):
                self._archive.close()
                return True
            elif self._archive is not None:
                sys.stdout.write("Wrong file type: " + str(type(self._archive)) + '\n')
                return False
            else:
                return True

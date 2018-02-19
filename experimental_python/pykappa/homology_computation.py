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

# Note: The script has to be called from sage itself
# sage -python /path/to/script.py

import inspect
import itertools
import sys
import time

from radial_cell import *
from parallel_cell import *
from misc_usability_stuff import *

def compute_homology(g=1, m=2, parametrized='True', verbose='True', homchain_file=None):
    # Setup coefficients and other variables.
    if g < 0 or (m < 0 if parametrized else m <= 0) :
        sys.stdout.write(
            'The genus has to be non-negative and the number of punctures has to be positive.\n'
            'Got g = ' + str(g) + '; m = ' + str(m) + '.\n'
        )
        sys.stdout.write('ABROATING!\n')
        sys.stdout.flush()

        raise ValueError()

    verbose = True if verbose == 'True' else False
    parametrized = True if parametrized == 'True' else False
    open_file = None
    if homchain_file is not None:
        sys.stdout.write('Opening file \'' + homchain_file + '\' ... ')
        sys.stdout.flush()
        try:
            open_file = open(homchain_file, 'wb')
        except:
            sys.stdout.write(' ABROATING! There was an error while writing the homchain file ' + homchain_file + '.\n')
            frameinfo = inspect.getframeinfo(inspect.currentframe())
            e, p, t = sys.exc_info()
            sys.stdout.write(str(frameinfo.filename) + ' ' + str(frameinfo.lineno) + '\n' + str(e) + ' ' + str(p) + '\n')
            sys.stdout.flush()
        else:
            sys.stdout.write('Done.\n')
            sys.stdout.flush()

    # Setup other variables.
    next_basis = {}
    dict_chaincomplex = {}
    degree = 4 * g + 2 * (m - 1) if parametrized == False else 4 * g + 2 * m
    starting_time = None

    # Load top cells.
    if verbose:
        sys.stdout.write(
            'We construct the cellular complex for g = ' + str(g) + ' and m = ' + str(m) + '.\n' +
            'Then we save the chain complex in chomp representation to \'' + homchain_file + '\'.\n\n'
        )
        sys.stdout.write('Loading top cells ... ')
        sys.stdout.flush()
        starting_time = time.clock()
    if g > 0 or m > 0:
        LoadTopCells = LoadTopRadialCellTau(g, m) if parametrized == False else LoadTopParallelCellTau(g, m)
        with LoadTopCells as rho_archive:
            # Check wether the archive exists or not.
            if rho_archive is None:
                sys.stdout.write("Loading top cells failed. We assume you did not create run our program 'create_and_store_top_cells.py' that creates and stores top cells for a given number hh = 2g+m-1.\n")
                sys.stdout.write('\n\n\n      ABROATING!\n\n\n')
                sys.stdout.flush()
                raise RuntimeError()

            for rho in rho_archive:
                cell = RadialCell(2*g+m-1, rho) if parametrized == False else ParallelCell(2*g+m, rho)
                next_basis[cell] = next_basis.get(cell, len(next_basis))
    if verbose:
        sys.stdout.write('Done. Duration = ' + str(time.clock() - starting_time) + '\n')

    # Compute the first differential.
    if verbose:
        sys.stdout.write('Computing the differential D_' + str(degree) + ' ... ')
        sys.stdout.flush()
        starting_time = time.clock()

    next_basis, num_rows, num_cols, bdry_matrix_dict = compute_faces_matrix(next_basis)
    dict_chaincomplex[degree] = (num_rows, num_cols, bdry_matrix_dict)

    degree -= 1

    if verbose:
        sys.stdout.write('Done. Duration = ' + str(time.clock() - starting_time) + '\n')

    # Compute the remaining differentials.
    while (next_basis is not None) and (len(next_basis) > 0):
        if verbose:
            sys.stdout.write('Computing the differential D_' + str(degree) + ' ... ')
            sys.stdout.flush()
            starting_time = time.clock()

        next_basis, num_rows, num_cols, bdry_matrix_dict = compute_faces_matrix(next_basis)
        dict_chaincomplex[degree] = (num_rows, num_cols, bdry_matrix_dict)

        if verbose:
            sys.stdout.write('Done. Duration = ' + str(time.clock() - starting_time) + '\n')
        degree -= 1

    # Construct the chain complex from the dictionaries.
    if verbose:
        sys.stdout.write('Constructing the chain complex ... ')
        sys.stdout.flush()
        starting_time = time.clock()

    sys.stdout.write('Writing chain complex representation to file \'' + homchain_file + '\' : \n')
    sys.stdout.flush()
    try:
        write_chaincomplex_to_chomp_representation(open_file, dict_chaincomplex, verbose)
    except:
        sys.stdout.write(' ABROATING! There was an error while writing the homchain file ' + homchain_file + '.\n')
        frameinfo = inspect.getframeinfo(inspect.currentframe())
        e, p, t = sys.exc_info()
        sys.stdout.write(str(frameinfo.filename) + ' ' + str(frameinfo.lineno) + '\n' + str(e) + ' ' + str(p) + '\n')
        sys.stdout.flush()
        try:
            open_file.close()
        except:
            pass
        raise
    else:
        try:
            open_file.close()
        except:
            raise
        return

def compute_faces_matrix(cells):
    # Setup variables.
    next_basis = {}
    matrix_dict = {}
    num_rows = 0
    num_cols = 0

    # Return empty stuff if input is empty.
    if cells is None or len(cells) == 0:
        return next_basis, num_rows, num_cols, matrix_dict

    # Get the degree.
    # The 'first' element of a dictionary is given by cells.iterkeys().next()
    degree = cells.iterkeys().next().degree()
    h = cells.iterkeys().next().get_h()
    if degree == 0:
        num_cols = len(cells)
        return next_basis, num_rows, num_cols, matrix_dict

    # Compute boundaries cell by cell.
    for cell in cells:
        cell_idx = cells[cell]
        # compute the boundary of the given cell
        coefficients = {}

        for s in itertools.product( *[ [i for i in range(1, q+1, 1)] for q in range(1, h+1, 1) ] ):
            current_basis = cell.get_clean_copy()
            norm_preserved = True
            parity = ((h * (h + 1)) / 2) % 2

            for q, s_q in enumerate(s):
                parity += s_q
                if current_basis.phi(q+1, s_q) is False:
                    norm_preserved = False
                    break
            parity = 1 if parity % 2 == 0 else -1

            if norm_preserved:
                sign = 1
                or_sign = current_basis.orientation_sign()

                for i in range(1, degree, 1):
                    face = current_basis.get_clean_copy()
                    if face.d_hor(i) == True:
                        face_idx = next_basis[face] = next_basis.get(face, len(next_basis))
                        coefficients[face_idx] = coefficients.get(face_idx, 0) + parity*sign*or_sign[i]
                    sign = -sign

        # store the column in the dictionary
        for face_idx, coeff in coefficients.items():
            matrix_dict[face_idx, cell_idx] = coeff


    # Compute number of columns and rows.
    num_cols = len(cells)
    num_rows = len(next_basis)

    # Done.
    return next_basis, num_rows, num_cols, matrix_dict


def write_chaincomplex_to_chomp_representation(open_file, diffs, verbose):
    # This function is an alternation of the SageMath function
    # "def _chomp_repr_(self)" in "SageMath/local/lib/python2.7/site-packages/sage/homology/chain_complex.py"
    # This provides an enormous speed up.

    if len(diffs) == 0:
        diffs = {0: (0, 0, [])}
    maxdim = max(diffs)
    mindim = 0
    s = "chain complex\n\nmax dimension = %s\n\n" % (maxdim - mindim)
    try:
        open_file.write(s)
    except:
        raise

    for i in range(0, maxdim + 1):
        if verbose:
            sys.stdout.write("  Dimension " + str(i) + " ... ")
            sys.stdout.flush()
            starting_time = time.clock()

        s = "dimension %s\n" % i
        num_rows, num_cols, mat = diffs.get(i, (0, 0, []))

        # We sort the matrix dictionary 'diffs' of the form
        #     (row, col) -> coeff
        # by the columns.
        # This is done by sort diff by the value of the second entry of the tuple (row, col) which is given
        # by the lambda function that gives the second value.
        sorted_keys = sorted(mat, key=lambda s: s[1])
        current_column_index = None
        idx = None
        non_zero_column = None
        # Print matrix if the matrix is not zero
        if len(sorted_keys) > 0:
            current_column_index = -1
            idx = -1
            for key in sorted_keys:
                # Iterate through columns.
                # Check if column is changed.
                if current_column_index != key[1]:
                    current_column_index = key[1]
                    idx += 1
                    # Finalize last column. Observe that we use True, False check in order to not finalize the very first differential.
                    if non_zero_column is True:
                        s += "\n"
                    elif non_zero_column is False:
                        s += "0\n"
                    s += "   boundary a%s = " % (idx + 1)
                    non_zero_column = False

                coeff = mat[key]
                if coeff > 0:
                    s += "+ %s * a%s " % (coeff, key[0] + 1)
                    non_zero_column = True
                elif coeff < 0:
                    s += "- %s * a%s " % (-coeff, key[0] + 1)
                    non_zero_column = True

            if non_zero_column == False:
                s += "0"
        else:
            # Print zero matrix if necessary.
            for idx in range(num_cols):
                s += "   boundary a%s = 0\n" % (idx + 1)
        s += "\n"

        try:
            open_file.write(s)
        except:
            raise

        if verbose:
            sys.stdout.write('Done. Duration = ' + str(time.clock() - starting_time) + '\n')
            sys.stdout.flush()


def main(g=0, m=7, parametrized=None, result_file=None, homchain_file=None):
    # Tee the output.
    tee = Tee(result_file, 'a')

    compute_homology(g, m, parametrized, 'True', homchain_file)
    sys.stdout.write('The Homology of the complex should be computed using the program \'homchain_gmp\'.\n')
    sys.stdout.flush()


main(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5])

#!/usr/bin/env python

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

import argparse
import cPickle as pickle
import gzip
import inspect
import sys
import time

import pykappa

# Store all top cells for which h = 2g+m-1.
def store_top_cells(h, parametrized):
    if parametrized:
        raise ValueError('Parametrization is not yet working.')

    if h < 1:
        print("h = 2g+m-1 is to small.")
        return

    sys.stdout.write('Creating and storing cells for 2g+m-1=' + str(h) + ' ... ' + '\r')
    sys.stdout.flush()
    starting_time = time.clock()

    # m and archive as a function of g
    m_of = []
    archive_of = []
    num_cells_of = ((h + 2) / 2)*[0]
    for g in range((h + 2) / 2):
        # exclude m = 0
        if h - 2 * g + 1 <= 0:
            continue
        m_of.append(h - 2 * g + 1)

        # Create directories.
        try:
            pykappa.create_directories()
            archive_of.append(gzip.GzipFile('./data/top_cell_g_' + str(g) + '_m_' + str(m_of[g]) + '.bz2', 'wb'))
        except:
            # Print the error.
            frameinfo = inspect.getframeinfo(inspect.currentframe())
            print(frameinfo.filename, frameinfo.lineno)
            e, p, t = sys.exc_info()
            print(e, p)
            return

    # iterate through all pairs
    def list_of_pairs(lst):
        if len(lst) < 2:
            yield lst
            return
        for i in range(1, len(lst)):
            pair = (lst[0], lst[i])
            for remaining in list_of_pairs(lst[1:i] + lst[i + 1:]):
                yield remaining + [pair]

    taus = list(list_of_pairs(range(2 * h, 0, -1)))
    l = len(taus)
    for n, tau in enumerate(taus):
        if n % 5000 == 0:
            sys.stdout.write(
                'Creating and storing cells for 2g+m-1=' + str(h) + ' ... ' + "{:.0%}".format(float(n) / l) + '\r')
            sys.stdout.flush()

        # Create Cell from permutation. This cell is a top cell by construction.
        cell = pykappa.RadialCell(h, tau)
        g = cell.get_g()
        pickle.dump(tau, archive_of[g])
        num_cells_of[g] += 1

    sys.stdout.write('Creating and storing cells for 2g+m-1=' + str(h) + ' ... Done. Number of cells: ' + ''.join( ['g='+ str(g) + ' -> ' + str(num_cells_of[g]) + '. ' for g in range((h+2)/2)] ) + 'Duration = ' + str(time.clock() - starting_time) + '\n')
    sys.stdout.flush()


def main():
    # check for correct version.
    major, minor = sys.version_info[0], sys.version_info[1]
    if major < 2 or (major == 2 and minor < 7):
        raise "Python >= 2.7 is required for argument parsing."

    # Use the argparse library. The library optparse is deprecated since version 2.7.
    # Compare the documentation: https://docs.python.org/2/library/argparse.html

    # Create the argument parser.
    # Note: Config files can be processed. In order to do so, we have to give fromfile_prefix_chars='_' with _ a symbol.
    parser = argparse.ArgumentParser(
        add_help=True,
        fromfile_prefix_chars='@',
        description='Caches the top cells.'
    )

    # Provide all arguments.
    # Note: we provide nargs=N and the N arguments from the command line will be gathered together into a list.
    # Thus, we supress nargs.
    parser.add_argument('-hh', required=True,   action='store', type=int, dest='h', metavar='arg', help='All top cells with h = 2g+m-1 = 2*genus + number of outgoing boundaries - 1.')
    parser.add_argument('-p', '--parametrized', action='store_true',      dest='parametrized',     help='Consider parametrized boundary.', default=False)
    args = vars(parser.parse_args())

    # Write Preamble.
    pre, valid = pykappa.preamble()
    sys.stdout.write(pre + '\n')
    sys.stdout.flush()

    # Stop of something went wrong.
    if valid == False:
        print("Could not initialize everything. Abroating.")
        return 1

    # Start actual computation.
    for h in range(1, args['h'] + 1):
        store_top_cells(h, args['parametrized'])

    # Cleanup screen with new lines. This is not written to the result file.
    sys.stdout.write('\n\n\n')
    sys.stdout.flush()


if __name__ == "__main__":
    main()

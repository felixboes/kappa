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
import subprocess
import sys

import pykappa

def call(g=1, m=2, parametrized=None, result_file=None, homcain_file=None):
    script_path = './pykappa/homology_computation.py'
    sys.stdout.write(
        'Calling python ' + script_path + ' ' + str(g) + ' ' + str(m) + ' ' + str(parametrized) + ' ' + str(result_file) +
        ' ' + str(homcain_file[0])  + '\n')
    sys.stdout.flush()
    cmd = ['python', script_path, str(g), str(m), str(parametrized), str(result_file), homcain_file[0]]
    subprocess.call(cmd)

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
        description='Compute the homology of the compactification of the unilevel radial slit domains aka sulivan diagrams.'
    )

    # Provide all arguments.
    # Note: we provide nargs=N and the N arguments from the command line will be gathered together into a list.
    # Thus, we supress nargs.
    parser.add_argument('-g', '--gen', required=True, action='store', type=int, dest='g', metavar='arg', help='The genus of the Riemann surfaces')
    parser.add_argument('-m', '--pun', required=True, action='store', type=int, dest='m', metavar='arg', help='The number of punctures of the Riemann surfaces')
    parser.add_argument('-p',                         action='store_true',      dest='parametrized',     help='If set, we use the parallel model', default=False)
    parser.add_argument('-s',          required=True, action='store', nargs=1,  dest='homcain_file', metavar='file', help='Store homchain file.')

    args = vars(parser.parse_args())

    # The name of the results file.
    args['result_file'] = './results/' + ''.join( [str(param).replace(' ', '_').replace('/', '_') for param in sys.argv if str(param)]).strip('.').lstrip('_')

    # Create directories.
    pykappa.create_directories()

    # Tee the output to the result file.
    tee = pykappa.Tee(args['result_file'], 'w')

    # Write Preamble.
    pre, valid = pykappa.preamble()
    sys.stdout.write(pre + '\n')
    sys.stdout.flush()

    # Since Sage uses file descriptors, we cannot use the Tee class directly.
    # That's unfortunate...
    tee.stop_logging()

    # Stop of something went wrong.
    if valid == False:
        print
        "Could not initialize everything. Abroating."
        return 1

    # Start actual computation.
    call(**args)

    # Cleanup screen with new lines. This is not written to the result file.
    sys.stdout.write('\n\n\n')
    sys.stdout.flush()


if __name__ == "__main__":
    main()

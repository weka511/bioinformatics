#!/usr/bin/env python

#   Copyright (C) 2023 Simon Crase

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

'''Generate shell script to execute all tests'''

from argparse import ArgumentParser
from glob     import glob
from os.path  import basename


def found_test(file_name):
    '''
    found_test

    Verify that file contains unit tests
    '''
    try:
        with open(file_name) as file:
            for line in file:
                if line.find('TestCase')>-1:
                    return True
    except UnicodeDecodeError:
        pass

    return False

def mark_executable(out):
    '''
    mark_executable
    Mark output file executable -- see answer to
    https://superuser.com/questions/778313/cygwins-chmod-behaves-as-working-but-it-does-not-work
    '''
    out.write('#!\n')

if __name__=='__main__':
    parser = ArgumentParser(__doc__)
    parser.add_argument('--exclude', nargs = '*', default=[],             help='List of files to be excluded')
    parser.add_argument('--out',     nargs = 1,   default='exectests.py', help='Name of batch file to be output')
    args    = parser.parse_args()
    exclude = args.exclude+[basename(__file__)] # Exclude this file as well as the ones the user specified
    with open(args.out,'w') as out:
        mark_executable(out)
        for file_name in glob('*.py'):
            if file_name not in exclude and found_test(file_name):
                out.write(f'echo {file_name}\n')
                out.write(f'./{file_name}\n')


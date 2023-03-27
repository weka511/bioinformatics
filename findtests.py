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

'''Generates shell script to execute all tests'''

from argparse import ArgumentParser
from glob     import glob
from os.path import basename
from time import time
from helpers import read_strings

def found_test(file_name):
    try:
        with open(file_name) as file:
            for line in file:
                if line.find('TestCase')>-1:
                    return True
    except UnicodeDecodeError:
        pass
    return False

if __name__=='__main__':
    moi = basename(__file__)
    with open('tests.sh','w') as out:
        for file_name in glob('*.py'):
            if moi != file_name and found_test(file_name):
                out.write(f'echo {file_name}\n')
                out.write(f'./{file_name}\n')


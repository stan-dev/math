#!/usr/bin/python3
#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from subprocess import run
import argparse
from typing import List
from os import path

# Helper to run the batch inserts example

_this_dir = path.abspath(path.dirname(path.realpath(__file__)))

class _Runner:
    def __init__(self, exe: str, host: str) -> None:
        self._exe = exe
        self._host = host
    
    def run(self, fname: str) -> None:
        json_path = path.join(_this_dir, fname)
        cmdline = [self._exe, 'example_user', 'example_password', self._host, json_path]
        print(' + ', ' '.join(cmdline))
        run(cmdline, check=True)


def main():
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('executable')
    parser.add_argument('host')
    args = parser.parse_args()

    # Build a runner
    runner = _Runner(args.executable, args.host)

    # Run the example with several combinations
    runner.run('employees_single.json')
    runner.run('employees_multiple.json')


if __name__ == '__main__':
    main()

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

# Helper to run the batch updates example

class _Runner:
    def __init__(self, exe: str, host: str) -> None:
        self._exe = exe
        self._host = host
    
    def run(self, updates: List[str]) -> None:
        employee_id = '1' # Guaranteed to exist
        cmdline = [self._exe, 'example_user', 'example_password', self._host, employee_id] + updates
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
    runner.run(['--salary=40000'])
    runner.run(['--company-id=HGS', '--last-name=Alice'])
    runner.run(['--salary=25000', '--company-id=AWC', '--first-name=John', '--last-name=Doe'])


if __name__ == '__main__':
    main()

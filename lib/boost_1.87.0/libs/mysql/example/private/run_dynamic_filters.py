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

# Helper to run the dynamic filters example

class _Runner:
    def __init__(self, exe: str, host: str) -> None:
        self._exe = exe
        self._host = host
    
    def run(self, opts: List[str]) -> None:
        cmdline = [self._exe, 'example_user', 'example_password', self._host] + opts
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
    runner.run(['--company-id=HGS'])
    runner.run(['--company-id=AWC', '--last-name=Alice'])
    runner.run(['--min-salary=25000', '--first-name=Bob', '--order-by=salary'])
    runner.run(['--company-id=AWC', '--first-name=Underpaid', '--last-name=Intern', '--min-salary=1', '--order-by=salary'])

if __name__ == '__main__':
    main()

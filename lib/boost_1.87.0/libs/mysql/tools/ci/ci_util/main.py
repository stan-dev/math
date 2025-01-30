#!/usr/bin/python3
#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from pathlib import Path
from typing import Union
import os
import argparse
from .common import IS_WINDOWS, BOOST_ROOT
from .cmake import cmake_build, cmake_noopenssl_build, cmake_nointeg_build, find_package_b2_test
from .b2 import b2_build
from .docs import docs_build


def _add_to_path(path: Path) -> None:
    sep = ';' if IS_WINDOWS  else ':'
    os.environ['PATH'] = '{}{}{}'.format(path, sep, os.environ["PATH"])


def _str2bool(v: Union[bool, str]) -> bool:
    if isinstance(v, bool):
        return v
    elif v == '1':
        return True
    elif v == '0':
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def _deduce_boost_branch() -> str:
    # Are we in GitHub Actions?
    if os.environ.get('GITHUB_ACTIONS') is not None:
        ci = 'GitHub Actions'
        ref = os.environ.get('GITHUB_BASE_REF', '') or os.environ.get('GITHUB_REF', '')
        res = 'master' if ref == 'master' or ref.endswith('/master') else 'develop'
    elif os.environ.get('DRONE') is not None:
        ref = os.environ.get('DRONE_BRANCH', '')
        ci = 'Drone'
        res = 'master' if ref == 'master' else 'develop'
    else:
        ci = 'Unknown'
        ref = ''
        res = 'develop'
    
    print('+  Found CI {}, ref={}, deduced branch {}'.format(ci, ref, res))

    return res


# Adds any DB args for builds requiring these
def _add_db_args(subp: argparse.ArgumentParser) -> None:
    subp.add_argument('--server-host', default='127.0.0.1')


# Fuzzing uses some Python features only available in newer CIs,
# so we use a dynamic import to avoid failures in older CIs
def _do_fuzz_build(**kwargs):
    from .fuzz import fuzz_build
    fuzz_build(**kwargs)


def main():
    # Main parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-dir', type=Path, required=True)
    parser.add_argument('--boost-branch', default=None) # None means "let this script deduce it"
    subparsers = parser.add_subparsers()

    # b2
    subp = subparsers.add_parser('b2', help='B2 build')
    _add_db_args(subp)
    subp.add_argument('--toolset', default='clang')
    subp.add_argument('--cxxstd', default='20')
    subp.add_argument('--variant', default='release')
    subp.add_argument('--stdlib', choices=['native', 'libc++'], default='native')
    subp.add_argument('--address-model', choices=['32', '64'], default='64')
    subp.add_argument('--separate-compilation', type=_str2bool, default=True)
    subp.add_argument('--use-ts-executor', type=_str2bool, default=False)
    subp.add_argument('--disable-local-sockets', default=None)
    subp.add_argument('--address-sanitizer', type=_str2bool, default=False)
    subp.add_argument('--undefined-sanitizer', type=_str2bool, default=False)
    subp.add_argument('--coverage', type=_str2bool, default=False)
    subp.add_argument('--valgrind', type=_str2bool, default=False)
    subp.add_argument('--fail-if-no-openssl', type=_str2bool, default=True)
    subp.set_defaults(func=b2_build)

    # cmake
    subp = subparsers.add_parser('cmake', help='CMake build')
    _add_db_args(subp)
    subp.add_argument('--generator', default='Ninja')
    subp.add_argument('--cmake-build-type', choices=['Debug', 'Release', 'MinSizeRel'], default='Debug', dest='build_type')
    subp.add_argument('--build-shared-libs', type=_str2bool, default=True)
    subp.add_argument('--cxxstd', default='20')
    subp.add_argument('--install-test', type=_str2bool, default=True)
    subp.set_defaults(func=cmake_build)

    # cmake without openssl
    subp = subparsers.add_parser('cmake-noopenssl', help='CMake build without OpenSSL')
    subp.add_argument('--generator', default='Ninja')
    subp.set_defaults(func=cmake_noopenssl_build)

    # cmake without integratin tests
    subp = subparsers.add_parser('cmake-nointeg', help='CMake build without integration tests')
    subp.add_argument('--generator', default='Ninja')
    subp.set_defaults(func=cmake_nointeg_build)

    # find_package with b2 distribution
    subp = subparsers.add_parser('find-package-b2', help='find_package with b2 distribution test')
    subp.add_argument('--generator', default='Ninja')
    subp.set_defaults(func=find_package_b2_test)

    # fuzz
    subp = subparsers.add_parser('fuzz', help='Fuzzing')
    _add_db_args(subp)
    subp.set_defaults(func=_do_fuzz_build)

    # docs
    subp = subparsers.add_parser('docs', help='Docs build')
    subp.set_defaults(func=docs_build)

    # Parse the arguments
    args = parser.parse_args()

    # Adjust Boost branch
    if args.boost_branch is None:
        args.boost_branch = _deduce_boost_branch()
    
    # Common code
    _add_to_path(BOOST_ROOT)
    
    # Call the build function, removing the func attribute
    args.func(**{name: value for name, value in vars(args).items() if name != 'func'})

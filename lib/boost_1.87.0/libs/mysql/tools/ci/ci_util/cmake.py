#!/usr/bin/python3
#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from pathlib import Path
import os
from typing import Optional, Dict
from .common import run, BOOST_ROOT, IS_WINDOWS, mkdir_and_cd
from .db_setup import db_setup
from .install_boost import install_boost


def _cmake_prefix_path(extra_path: Optional[Path] = None) -> str:
    prefix_path_args = [extra_path, 'C:\\openssl-64' if IS_WINDOWS else None]
    return ';'.join(str(p) for p in prefix_path_args if p is not None)


def _cmake_bool(value: bool) -> str:
    return 'ON' if value else 'OFF'


class _CMakeRunner:
    def __init__(self, generator: str, build_type: str) -> None:
        self._generator = generator
        self._build_type = build_type
        os.environ['CMAKE_BUILD_PARALLEL_LEVEL'] = '4'


    def configure(self, source_dir: Path, binary_dir: Path, variables: Dict[str, str]) -> None:
        mkdir_and_cd(binary_dir)
        run(
            [
                'cmake',
                '-G',
                self._generator,
                '-DCMAKE_BUILD_TYPE={}'.format(self._build_type),
            ] +
            ['-D{}={}'.format(name, value) for name, value in variables.items()] +
            [str(source_dir)]
        )


    def build(self, target: str) -> None:
        run(['cmake', '--build', '.', '--target', target, '--config', self._build_type])
    

    def build_all(self) -> None:
        run(['cmake', '--build', '.', '--config', self._build_type])


    def ctest(self, no_tests_error: bool = True) -> None:
        run(
            ['ctest', '--output-on-failure'] +
            (['--no-tests=error'] if no_tests_error else []) +
            ['--build-config', self._build_type]
        )


# Regular CMake build
def cmake_build(
    source_dir: Path,
    boost_branch: str,
    generator: str,
    build_shared_libs: bool,
    build_type: str,
    cxxstd: str,
    server_host: str,
    install_test: bool,
) -> None:
    # Config
    cmake_distro = Path(os.path.expanduser('~')).joinpath('cmake-distro')
    test_folder = BOOST_ROOT.joinpath('libs', 'mysql', 'test', 'cmake_test')
    runner = _CMakeRunner(generator=generator, build_type=build_type)

    # Get Boost
    install_boost(
        source_dir=source_dir,
        boost_branch=boost_branch,
    )

    # Setup DB
    db_setup(source_dir, server_host)

    # Build the library, run the tests, and install, as the Boost superproject does
    bin_dir = BOOST_ROOT.joinpath('__build')
    runner.configure(
        source_dir=BOOST_ROOT,
        binary_dir=bin_dir,
        variables={
            'CMAKE_PREFIX_PATH': _cmake_prefix_path(),
            'BOOST_INCLUDE_LIBRARIES': 'mysql',
            'BUILD_SHARED_LIBS': _cmake_bool(build_shared_libs),
            'CMAKE_INSTALL_PREFIX': str(cmake_distro),
            'BUILD_TESTING': 'ON',
            'CMAKE_INSTALL_MESSAGE': 'NEVER',
            'BOOST_MYSQL_INTEGRATION_TESTS': 'ON',
            **({ 'CMAKE_CXX_STANDARD': cxxstd } if cxxstd else {})
        }
    )
    runner.build(target='tests')
    runner.ctest()
    runner.build(target='install')

    # The library can be consumed using add_subdirectory
    runner.configure(
        source_dir=test_folder,
        binary_dir=test_folder.joinpath('__build_add_subdirectory'),
        variables={
            'CMAKE_PREFIX_PATH': _cmake_prefix_path(),
            'BOOST_CI_INSTALL_TEST': 'OFF',
            'BUILD_SHARED_LIBS': _cmake_bool(build_shared_libs)
        }
    )
    runner.build_all()
    runner.ctest()

    # The library can be consumed using find_package on a Boost distro built by cmake.
    # Generating a modular installation requires CMake 3.13+, so this test can be disabled
    # for jobs running old cmake versions
    if install_test:
        runner.configure(
            source_dir=test_folder,
            binary_dir=test_folder.joinpath('__build_find_package'),
            variables={
                'BOOST_CI_INSTALL_TEST': 'ON',
                'BUILD_SHARED_LIBS': _cmake_bool(build_shared_libs),
                'CMAKE_PREFIX_PATH': _cmake_prefix_path(cmake_distro)
            }
        )
        runner.build_all()
        runner.ctest()


# Check that we bail out correctly when no OpenSSL is available
def cmake_noopenssl_build(
    source_dir: Path,
    boost_branch: str,
    generator: str
):
    # Config
    runner = _CMakeRunner(generator=generator, build_type='Release')

    # Get Boost
    install_boost(
        source_dir=source_dir,
        boost_branch=boost_branch,
    )

    # Build Boost and install. This won't build the library if there is no OpenSSL,
    # but shouldn't cause any error.
    runner.configure(
        source_dir=BOOST_ROOT,
        binary_dir=BOOST_ROOT.joinpath('__build'),
        variables={
            'BOOST_INCLUDE_LIBRARIES': 'mysql',
            'BOOST_MYSQL_INTEGRATION_TESTS': 'ON',
            'BUILD_TESTING': 'ON',
            'CMAKE_INSTALL_MESSAGE': 'NEVER',
        }
    )
    runner.build(target='tests')
    runner.ctest(no_tests_error=False)
    runner.build(target='install')


# Check that disabling integration tests works
def cmake_nointeg_build(
    source_dir: Path,
    boost_branch: str,
    generator: str
):
    # Config
    runner = _CMakeRunner(generator=generator, build_type='Release')

    # Get Boost
    install_boost(
        source_dir=source_dir,
        boost_branch=boost_branch,
    )

    # Build the library and run the tests, as the Boost superproject does
    bin_dir = BOOST_ROOT.joinpath('__build')
    runner.configure(
        source_dir=BOOST_ROOT,
        binary_dir=bin_dir,
        variables={
            'CMAKE_PREFIX_PATH': _cmake_prefix_path(),
            'BOOST_INCLUDE_LIBRARIES': 'mysql',
            'BUILD_TESTING': 'ON'
        }
    )
    runner.build(target='tests')
    runner.ctest()
    runner.build(target='install')


# Check that CMake can consume a Boost distribution created by b2
def find_package_b2_test(
    source_dir: Path,
    boost_branch: str,
    generator: str
) -> None:
    # Config
    b2_distro = Path(os.path.expanduser('~')).joinpath('b2-distro')
    prefix_path = _cmake_prefix_path(b2_distro)
    runner = _CMakeRunner(generator=generator, build_type='Release')

    # Get Boost
    install_boost(
        source_dir=source_dir,
        boost_branch=boost_branch,
    )

    # Generate a b2 Boost distribution
    run([
        'b2',
        '--prefix={}'.format(b2_distro),
        '--with-charconv',
        '-d0',
        'install'
    ])

    # Check that the library can be consumed using find_package on the distro above
    test_dir = BOOST_ROOT.joinpath('libs', 'mysql', 'test', 'cmake_b2_test')
    runner.configure(
        source_dir=test_dir,
        binary_dir=test_dir.joinpath('__build'),
        variables={
            'CMAKE_PREFIX_PATH': prefix_path,
            'BUILD_TESTING': 'ON'
        }
    )
    runner.build_all()
    runner.ctest()

    # Same as the above, but for separate-build mode
    test_dir = BOOST_ROOT.joinpath('libs', 'mysql', 'test', 'cmake_b2_separate_compilation_test')
    runner.configure(
        source_dir=test_dir,
        binary_dir=test_dir.joinpath('__build'),
        variables={
            'CMAKE_PREFIX_PATH': prefix_path,
            'BUILD_TESTING': 'ON'
        }
    )
    runner.build_all()
    runner.ctest()


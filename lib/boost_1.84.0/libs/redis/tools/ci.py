#!/usr/bin/python3

# Contains commands that are invoked by the CI scripts.
# Having this as a Python file makes it platform-independent.

from pathlib import Path
from typing import List, Union
import subprocess
import os
import stat
from shutil import rmtree, copytree, ignore_patterns
import argparse


# Variables
_is_windows = os.name == 'nt'
_home = Path(os.path.expanduser('~'))
_boost_root = _home.joinpath('boost-root')
_b2_distro = _home.joinpath('boost-b2-distro')
_cmake_distro = _home.joinpath('boost-cmake-distro')
_b2_command = str(_boost_root.joinpath('b2'))


# Utilities
def _run(args: List[str]) -> None:
    print('+ ', args, flush=True)
    subprocess.run(args, check=True)


def _mkdir_and_cd(path: Path) -> None:
    os.makedirs(str(path), exist_ok=True)
    os.chdir(str(path))


def _cmake_bool(value: bool) -> str:
    return 'ON' if value else 'OFF'


def _remove_readonly(func, path, _):
    os.chmod(path, stat.S_IWRITE)
    func(path)


# Parses a string into a boolean (for command-line parsing)
def _str2bool(v: Union[bool, str]) -> bool:
    if isinstance(v, bool):
        return v
    elif v == '1':
        return True
    elif v == '0':
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# Transforms a b2-like toolset into a compiler command suitable
# to be passed to CMAKE_CXX_COMPILER
def _compiler_from_toolset(toolset: str) -> str:
    if toolset.startswith('gcc'):
        return toolset.replace('gcc', 'g++')
    elif toolset.startswith('clang'):
        return toolset.replace('clang', 'clang++')
    elif toolset.startswith('msvc'):
        return 'cl'
    else:
        return toolset


# If we're on the master branch, we should use the Boost superproject master branch.
# Otherwise, use the superproject develop branch.
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


# Gets Boost (develop or master, depending on the CI branch we're operating on),
# with the required dependencies, and leaves it at _boost_root. Places our library,
# located under source_dir, under $BOOST_ROOT/libs. Also runs the bootstrap script so b2 is usable. 
def _setup_boost(
    source_dir: Path
) -> None:
    assert source_dir.is_absolute()
    assert not _boost_root.exists()
    lib_dir = _boost_root.joinpath('libs', 'redis')
    branch = _deduce_boost_branch()

    # Clone Boost
    _run(['git', 'clone', '-b', branch, '--depth', '1', 'https://github.com/boostorg/boost.git', str(_boost_root)])
    os.chdir(str(_boost_root))

    # Put our library inside boost root
    if lib_dir.exists():
        rmtree(str(lib_dir), onerror=_remove_readonly)
    copytree(
        str(source_dir),
        str(lib_dir),
        ignore=ignore_patterns('__build*__', '.git')
    )

    # Install Boost dependencies
    _run(["git", "config", "submodule.fetchJobs", "8"])
    _run(["git", "submodule", "update", "-q", "--init", "tools/boostdep"])
    _run(["python", "tools/boostdep/depinst/depinst.py", "--include", "example", "redis"])

    # Bootstrap
    if _is_windows:
        _run(['cmd', '/q', '/c', 'bootstrap.bat'])
    else:
        _run(['bash', 'bootstrap.sh'])
    _run([_b2_command, 'headers', '-d0'])


# Builds a Boost distribution using ./b2 install, and places it into _b2_distro.
# This emulates a regular Boost distribution, like the ones in releases
def _build_b2_distro(
    toolset: str
):
    os.chdir(str(_boost_root))
    _run([
        _b2_command,
        '--prefix={}'.format(_b2_distro),
        '--with-system',
        'toolset={}'.format(toolset),
        '-d0',
        'install'
    ])


# Builds a Boost distribution using cmake, and places it into _cmake_distro.
# It includes only our library and any dependency.
def _build_cmake_distro(
    generator: str,
    build_type: str,
    cxxstd: str,
    toolset: str,
    build_shared_libs: bool = False
):
    _mkdir_and_cd(_boost_root.joinpath('__build_cmake_test__'))
    _run([
        'cmake',
        '-G',
        generator,
        '-DBUILD_TESTING=ON',
        '-DCMAKE_CXX_COMPILER={}'.format(_compiler_from_toolset(toolset)),
        '-DCMAKE_BUILD_TYPE={}'.format(build_type),
        '-DCMAKE_CXX_STANDARD={}'.format(cxxstd),
        '-DBOOST_INCLUDE_LIBRARIES=redis',
        '-DBUILD_SHARED_LIBS={}'.format(_cmake_bool(build_shared_libs)),
        '-DCMAKE_INSTALL_PREFIX={}'.format(_cmake_distro),
        '-DBUILD_TESTING=ON',
        '-DBoost_VERBOSE=ON',
        '-DCMAKE_INSTALL_MESSAGE=NEVER',
        '..'
    ])
    _run(['cmake', '--build', '.', '--target', 'tests', '--config', build_type])
    _run(['ctest', '--output-on-failure', '--build-config', build_type])
    _run(['cmake', '--build', '.', '--target', 'install', '--config', build_type])


# Builds our CMake tests as a standalone project
# (BOOST_REDIS_MAIN_PROJECT is ON) and we find_package Boost.
# This ensures that all our test suite is run.
def _build_cmake_standalone_tests(
    generator: str,
    build_type: str,
    cxxstd: str,
    toolset: str,
    build_shared_libs: bool = False
):
    _mkdir_and_cd(_boost_root.joinpath('libs', 'redis', '__build_standalone__'))
    _run([
        'cmake',
        '-DBUILD_TESTING=ON',
        '-DCMAKE_CXX_COMPILER={}'.format(_compiler_from_toolset(toolset)),
        '-DCMAKE_PREFIX_PATH={}'.format(_b2_distro),
        '-DCMAKE_BUILD_TYPE={}'.format(build_type),
        '-DBUILD_SHARED_LIBS={}'.format(_cmake_bool(build_shared_libs)),
        '-DCMAKE_CXX_STANDARD={}'.format(cxxstd),
        '-G',
        generator,
        '..'
    ])
    _run(['cmake', '--build', '.'])


# Runs the tests built in the previous step
def _run_cmake_standalone_tests(
    build_type: str
):
    os.chdir(str(_boost_root.joinpath('libs', 'redis', '__build_standalone__')))
    _run(['ctest', '--output-on-failure', '--build-config', build_type, '--no-tests=error'])


# Tests that the library can be consumed using add_subdirectory()
def _run_cmake_add_subdirectory_tests(
    generator: str,
    build_type: str,
    cxxstd: str,
    toolset: str,
    build_shared_libs: bool = False
):
    test_folder = _boost_root.joinpath('libs', 'redis', 'test', 'cmake_subdir_test', '__build')
    _mkdir_and_cd(test_folder)
    _run([
        'cmake',
        '-G',
        generator,
        '-DCMAKE_CXX_COMPILER={}'.format(_compiler_from_toolset(toolset)),
        '-DBUILD_TESTING=ON',
        '-DCMAKE_BUILD_TYPE={}'.format(build_type),
        '-DBUILD_SHARED_LIBS={}'.format(_cmake_bool(build_shared_libs)),
        '-DCMAKE_CXX_STANDARD={}'.format(cxxstd),
        '..'
    ])
    _run(['cmake', '--build', '.', '--config', build_type])
    _run(['ctest', '--output-on-failure', '--build-config', build_type, '--no-tests=error'])


# Tests that the library can be consumed using find_package on a distro built by cmake
def _run_cmake_find_package_tests(
    generator: str,
    build_type: str,
    cxxstd: str,
    toolset: str,
    build_shared_libs: bool = False
):
    _mkdir_and_cd(_boost_root.joinpath('libs', 'redis', 'test', 'cmake_install_test', '__build'))
    _run([
        'cmake',
        '-G',
        generator,
        '-DCMAKE_CXX_COMPILER={}'.format(_compiler_from_toolset(toolset)),
        '-DBUILD_TESTING=ON',
        '-DCMAKE_BUILD_TYPE={}'.format(build_type),
        '-DBUILD_SHARED_LIBS={}'.format(_cmake_bool(build_shared_libs)),
        '-DCMAKE_CXX_STANDARD={}'.format(cxxstd),
        '-DCMAKE_PREFIX_PATH={}'.format(_cmake_distro),
        '..'
    ])
    _run(['cmake', '--build', '.', '--config', build_type])
    _run(['ctest', '--output-on-failure', '--build-config', build_type, '--no-tests=error'])


# Tests that the library can be consumed using find_package on a distro built by b2
def _run_cmake_b2_find_package_tests(
    generator: str,
    build_type: str,
    cxxstd: str,
    toolset: str,
    build_shared_libs: bool = False
):
    _mkdir_and_cd(_boost_root.joinpath('libs', 'redis', 'test', 'cmake_b2_test', '__build'))
    _run([
        'cmake',
        '-G',
        generator,
        '-DCMAKE_CXX_COMPILER={}'.format(_compiler_from_toolset(toolset)),
        '-DBUILD_TESTING=ON',
        '-DCMAKE_PREFIX_PATH={}'.format(_b2_distro),
        '-DCMAKE_BUILD_TYPE={}'.format(build_type),
        '-DBUILD_SHARED_LIBS={}'.format(_cmake_bool(build_shared_libs)),
        '-DCMAKE_CXX_STANDARD={}'.format(cxxstd),
        '-DBUILD_TESTING=ON',
        '..'
    ])
    _run(['cmake', '--build', '.', '--config', build_type])
    _run(['ctest', '--output-on-failure', '--build-config', build_type, '--no-tests=error'])


# Builds and runs the library tests using b2
def _run_b2_tests(
    variant: str,
    cxxstd: str,
    toolset: str
):
    os.chdir(str(_boost_root))
    _run([
        _b2_command,
        '--abbreviate-paths',
        'toolset={}'.format(toolset),
        'cxxstd={}'.format(cxxstd),
        'variant={}'.format(variant),
        '-j4',
        'libs/redis/test'
    ])


def main():
    # Command line parsing
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    subp = subparsers.add_parser('setup-boost')
    subp.add_argument('--source-dir', type=Path, required=True)
    subp.set_defaults(func=_setup_boost)

    subp = subparsers.add_parser('build-b2-distro')
    subp.add_argument('--toolset', default='gcc')
    subp.set_defaults(func=_build_b2_distro)

    subp = subparsers.add_parser('build-cmake-distro')
    subp.add_argument('--generator', default='Unix Makefiles')
    subp.add_argument('--build-type', default='Debug')
    subp.add_argument('--cxxstd', default='20')
    subp.add_argument('--toolset', default='gcc')
    subp.add_argument('--build-shared-libs', type=_str2bool, default=False)
    subp.set_defaults(func=_build_cmake_distro)

    subp = subparsers.add_parser('build-cmake-standalone-tests')
    subp.add_argument('--generator', default='Unix Makefiles')
    subp.add_argument('--build-type', default='Debug')
    subp.add_argument('--cxxstd', default='20')
    subp.add_argument('--toolset', default='gcc')
    subp.add_argument('--build-shared-libs', type=_str2bool, default=False)
    subp.set_defaults(func=_build_cmake_standalone_tests)

    subp = subparsers.add_parser('run-cmake-standalone-tests')
    subp.add_argument('--build-type', default='Debug')
    subp.set_defaults(func=_run_cmake_standalone_tests)

    subp = subparsers.add_parser('run-cmake-add-subdirectory-tests')
    subp.add_argument('--generator', default='Unix Makefiles')
    subp.add_argument('--build-type', default='Debug')
    subp.add_argument('--cxxstd', default='20')
    subp.add_argument('--toolset', default='gcc')
    subp.add_argument('--build-shared-libs', type=_str2bool, default=False)
    subp.set_defaults(func=_run_cmake_add_subdirectory_tests)

    subp = subparsers.add_parser('run-cmake-find-package-tests')
    subp.add_argument('--generator', default='Unix Makefiles')
    subp.add_argument('--build-type', default='Debug')
    subp.add_argument('--cxxstd', default='20')
    subp.add_argument('--toolset', default='gcc')
    subp.add_argument('--build-shared-libs', type=_str2bool, default=False)
    subp.set_defaults(func=_run_cmake_find_package_tests)

    subp = subparsers.add_parser('run-cmake-b2-find-package-tests')
    subp.add_argument('--generator', default='Unix Makefiles')
    subp.add_argument('--build-type', default='Debug')
    subp.add_argument('--cxxstd', default='20')
    subp.add_argument('--toolset', default='gcc')
    subp.add_argument('--build-shared-libs', type=_str2bool, default=False)
    subp.set_defaults(func=_run_cmake_b2_find_package_tests)

    subp = subparsers.add_parser('run-b2-tests')
    subp.add_argument('--variant', default='debug,release')
    subp.add_argument('--cxxstd', default='17,20')
    subp.add_argument('--toolset', default='gcc')
    subp.set_defaults(func=_run_b2_tests)

    # Actually parse the arguments
    args = parser.parse_args()

    # Invoke the relevant function (as defined by the func default), with
    # the command-line arguments the user passed us (we need to get rid
    # of the func property to match function signatures)
    # This approach is recommended by Python's argparse docs
    args.func(**{k: v for k, v in vars(args).items() if k != 'func'})


if __name__ == '__main__':
    main()

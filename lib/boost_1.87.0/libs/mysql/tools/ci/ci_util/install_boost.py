#!/usr/bin/python3
#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from shutil import rmtree, ignore_patterns, copytree
import stat
import sys
from pathlib import Path
import os
from .common import run, IS_WINDOWS, BOOST_ROOT

def _remove_readonly(func, path, _):
    os.chmod(path, stat.S_IWRITE)
    func(path)

def _copy_lib_to_boost(source_dir: Path):
    # Config
    supports_dir_exist_ok = sys.version_info.minor >= 8
    lib_dir = BOOST_ROOT.joinpath('libs', 'mysql')
    
    # Old versions of Python don't support dirs_exist_ok.
    # For these, we need to remove any old copies of our lib before copying
    if lib_dir.exists() and not supports_dir_exist_ok:
        rmtree(str(lib_dir), onerror=_remove_readonly)
    
    # Do the copying
    copytree(
        str(source_dir),
        str(lib_dir),
        ignore=ignore_patterns('__build*__', '.git'),
        **({ 'dirs_exist_ok': True } if supports_dir_exist_ok else {}) # type: ignore
    )


def install_boost(
    source_dir: Path,
    boost_branch: str,
    docs_install: bool = False
) -> None:
    assert source_dir.is_absolute()
    
    # If BOOST_ROOT already exists, this is a re-build.
    # Copy our library into libs/ and exit
    if BOOST_ROOT.exists():
        os.chdir(str(BOOST_ROOT))
        _copy_lib_to_boost(source_dir)
        return

    # Clone Boost
    run(['git', 'clone', '-b', boost_branch, '--depth', '1', 'https://github.com/boostorg/boost.git', str(BOOST_ROOT)])
    os.chdir(str(BOOST_ROOT))

    # Put our library inside boost root
    _copy_lib_to_boost(source_dir)

    # Install Boost dependencies
    submodules = [
        'libs/context',
        'tools/boostdep',
        'tools/boostbook',
        'tools/docca',
        'tools/quickbook'
    ] if docs_install else [
        'tools/boostdep'
    ]
    run(["git", "config", "submodule.fetchJobs", "8"])
    run(["git", "submodule", "update", "-q", "--init"] + submodules)

    if docs_install:
        run(['python', 'tools/boostdep/depinst/depinst.py', '../tools/quickbook'])
    else:
        run(["python", "tools/boostdep/depinst/depinst.py", "--include", "example", "mysql"])
    
    # Bootstrap
    if IS_WINDOWS:
        run(['cmd', '/q', '/c', 'bootstrap.bat'])
    else:
        run(['bash', 'bootstrap.sh'])
    run(['b2', 'headers'])

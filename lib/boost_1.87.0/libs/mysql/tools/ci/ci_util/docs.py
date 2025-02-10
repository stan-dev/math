#!/usr/bin/python3
#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from pathlib import Path
import os
from shutil import copytree, rmtree
from .common import BOOST_ROOT, run
from .install_boost import install_boost


def docs_build(
    source_dir: Path,
    boost_branch: str
):
    # Get Boost. This leaves us inside boost root
    install_boost(
        source_dir=source_dir,
        boost_branch=boost_branch,
        docs_install=True,
    )

    # Write the config file
    config_path = os.path.expanduser('~/user-config.jam')
    with open(config_path, 'wt') as f:
        f.writelines(['using doxygen ;\n', 'using boostbook ;\n'])

    # Run b2
    run(['b2', 'libs/mysql/doc//boostrelease'])

    # Copy the resulting docs into a well-known path
    output_dir = source_dir.joinpath('doc', 'html')
    if output_dir.exists():
        rmtree(output_dir)
    copytree(BOOST_ROOT.joinpath('libs', 'mysql', 'doc', 'html'), output_dir)

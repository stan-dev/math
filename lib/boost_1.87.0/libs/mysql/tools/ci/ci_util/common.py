#!/usr/bin/python3
#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from pathlib import Path
import os
from typing import List
import subprocess

REPO_BASE = Path(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', '..')).absolute()
BOOST_ROOT = Path(os.path.expanduser('~')).joinpath('boost-root')
IS_WINDOWS = os.name == 'nt'


def run(args: List[str]) -> None:
    print('+ ', args, flush=True)
    subprocess.run(args, check=True)


def mkdir_and_cd(path: Path) -> None:
    os.makedirs(str(path), exist_ok=True)
    os.chdir(str(path))

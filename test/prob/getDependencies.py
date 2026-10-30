#!/usr/bin/env python3
"""
This script determines which distributions may have changed
based on the files #included by their test. It then prints the distribution
fixtures (test/prob/<dist>/<fixture>_test.hpp) whose tests should run; feed
them to the top-level runTests.py script.

The dependencies come from the compiler with the flags of a configured CMake
build directory (runTests.py --configure-only), so they match the test build.
"""

from concurrent.futures import ThreadPoolExecutor
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
from typing import List, Set

ROOT = Path(__file__).resolve().parents[2]


def dependency_command(compile_commands: List[dict]) -> List[str]:
    """Turn the gtest_main compile command into a dependency listing command."""
    for entry in compile_commands:
        if entry["file"].endswith("gtest_main.cc"):
            arguments = entry.get("arguments") or shlex.split(entry["command"])
            break
    else:
        raise ValueError("No gtest_main.cc compile command; configure with tests enabled")
    command, skip = [], False
    for argument in arguments:
        if skip:
            skip = False
        elif argument in ("-o", "-MT", "-MF", "-MQ"):
            skip = True
        elif argument in ("-c", "-MD", "-MMD") or argument.endswith("gtest_main.cc"):
            continue
        else:
            command.append(argument)
    return command + ["-MM", "-MG", "-x", "c++"]


def parse_dependencies(output: str, root: Path) -> Set[str]:
    """Keep the stan/math headers from make-style dependency output."""
    words = output.replace("\\\n", " ").split(":", 1)[-1].split()
    deps = set()
    for word in words:
        path = Path(word) if Path(word).is_absolute() else root / word
        try:
            relative = Path(os.path.normpath(path)).relative_to(root).as_posix()
        except ValueError:
            continue
        if "stan/math" in relative and relative.endswith(".hpp"):
            deps.add(relative)
    return deps


def get_dependencies(command: List[str], fixture: str, root: Path) -> Set[str]:
    result = subprocess.run(command + [fixture], cwd=root, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"Listing dependencies of {fixture} failed:\n{result.stderr}")
    return parse_dependencies(result.stdout, root) | {fixture}


def get_changed(root: Path) -> Set[str]:
    return set(subprocess.run(["git", "diff", "--name-only", "origin/develop...HEAD"],
                              cwd=root, text=True, capture_output=True,
                              check=True).stdout.splitlines())


def main(arguments, root: Path) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pretend-all", action="store_true",
                        help="select every distribution fixture")
    parser.add_argument("--build-dir", default="build",
                        help="configured CMake build directory (default: build)")
    inputs = parser.parse_args(arguments)

    fixtures = sorted(path.relative_to(root).as_posix()
                      for path in (root / "test" / "prob").glob("*/*_test.hpp"))
    if inputs.pretend_all:
        selected = fixtures
    else:
        compile_commands = Path(inputs.build_dir) / "compile_commands.json"
        if not compile_commands.is_absolute():
            compile_commands = Path.cwd() / compile_commands
        command = dependency_command(json.loads(compile_commands.read_text()))
        changed = get_changed(root)
        with ThreadPoolExecutor(os.cpu_count()) as pool:
            all_deps = list(pool.map(lambda fixture: get_dependencies(command, fixture, root),
                                     fixtures))
        selected = []
        for fixture, deps in zip(fixtures, all_deps):
            intersection = changed & deps
            if intersection:
                print(f"{fixture} has {len(intersection)} changed dependencies out of "
                      f"{len(deps)} files #included.", file=sys.stderr)
                selected.append(fixture)
    if selected:
        print("\n".join(selected))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:], ROOT))

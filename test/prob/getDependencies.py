#!/usr/bin/env python3
"""
This script determines which distributions may have changed
based on the files #included by their test. It then prints a set of strings
which will run said tests when fed to the top-level runTests.py script
"""

from typing import Set
import sys
import subprocess
from pathlib import Path
import glob

def get_dependencies(file: Path) -> Set[str]:
    file_dot_d = file.with_suffix(".d")
    subprocess.run(["make", str(file_dot_d)], stdout=subprocess.DEVNULL)
    contents = file_dot_d.read_text()
    file_dot_d.unlink()

    deps = set(
        filter(
            lambda s: "stan/math" in s and s.endswith(".hpp"),
            map(str.strip, contents.split(" ")),
        )
    )
    return deps | {str(file)}


def get_changed() -> Set[str]:
    changed_files = subprocess.run(
        ["git", "diff", "--name-only", "origin/develop...HEAD"],
        text=True,
        capture_output=True,
    ).stdout.splitlines()

    return set(changed_files)


def add_tests_from_hpp(tests_to_run, test):
    pattern = str(test).replace("_test.hpp", "_0*_test.cpp")
    tests = glob.glob(pattern)
    tests_to_run.extend(tests)

if __name__ == "__main__":
    if Path("makefile") not in Path(".").iterdir():
        raise ValueError("getDependencies must be ran from the top-level repository")

    pretend = "--pretend-all" in sys.argv

    changed = get_changed()
    tests_to_run = []

    distribution_tests = Path("test", "prob")

    subprocess.run(["make", "generate-tests"], stdout=subprocess.DEVNULL)

    for dist in distribution_tests.iterdir():
        if not dist.is_dir():
            continue
        print(f"Considering {dist}.", file=sys.stderr)
        for test in dist.iterdir():
            if test.suffix != ".hpp": continue
            if pretend:
                add_tests_from_hpp(tests_to_run, test)
            else:
                deps = get_dependencies(test)
                intersection = changed & deps
                if intersection:
                    print(
                        f"{test} has {len(intersection)} changed dependencies out of {len(deps)} files #included.",
                        file=sys.stderr,
                    )
                    add_tests_from_hpp(tests_to_run, test)

    print("\n".join(set(tests_to_run)))

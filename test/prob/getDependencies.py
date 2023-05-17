#!/usr/bin/python
"""
This script determines which distributions may have changed
based on the files #included by their test. It then prints a set of strings
which will run said tests when fed to the top-level runTests.py script
"""

from typing import Set
import sys
import subprocess
from pathlib import Path


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


if __name__ == "__main__":
    if Path("makefile") not in Path(".").iterdir():
        raise ValueError("getDependencies must be ran from the top-level repository")

    pretend = "--pretend-all" in sys.argv

    changed = get_changed()
    tests_to_run = []

    distribution_tests = Path("test", "prob")

    for dist in distribution_tests.iterdir():
        if not dist.is_dir():
            continue
        for test in dist.iterdir():
            if pretend:
                tests_to_run.append(str(test).replace("_test.hpp", "_0*_test.cpp"))
            else:
                deps = get_dependencies(test)
                intersection = changed & deps
                if intersection:
                    print(
                        f"{test} has {len(intersection)} changed dependencies out of {len(deps)} files #included.",
                        file=sys.stderr,
                    )

                    tests_to_run.append(str(test).replace("_test.hpp", "_0*_test.cpp"))

    print("\n".join(tests_to_run))

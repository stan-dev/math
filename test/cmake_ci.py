#!/usr/bin/env python3
"""Select CI tests and generate expression sources for CMake builds."""

import argparse
from pathlib import Path
import re
import subprocess
import sys


def changed_files():
    return subprocess.run(
        ["git", "diff", "--name-only", "--diff-filter=d", "origin/develop...HEAD"],
        check=True, text=True, capture_output=True,
    ).stdout.splitlines()


def source_filter(paths):
    return "^(" + "|".join(re.escape(path) for path in sorted(set(paths))) + ")$" if paths else ""


def changed_units():
    # Preserve the existing runner's priority-test selection and size limit.
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from runTests import findChangedTests

    return source_filter([path + ".cpp" for path in findChangedTests(False)])


def changed_distributions(files, run_all=False):
    # Every distribution fixture includes the Stan umbrella headers. Library,
    # build-system, and shared test-header changes can affect all distributions.
    if run_all or any(
        path.startswith(("stan/", "lib/", "make/", "cmake/"))
        or path in ("CMakeLists.txt", "Jenkinsfile", "makefile")
        or (path.startswith("test/") and path.endswith(".hpp")
            and not re.fullmatch(r"test/prob/[^/]+/[^/]+_test\.hpp", path))
        or (path.startswith("test/prob/") and not path.endswith("_test.hpp"))
        for path in files
    ):
        return "^test/prob/"
    fixtures = [path for path in files
                if re.fullmatch(r"test/prob/[^/]+/[^/]+_test\.hpp", path)]
    return ("^(" + "|".join(re.escape(path[:-len("_test.hpp")])
                            for path in sorted(fixtures)) + ")_[0-9]+_generated_"
            if fixtures else "")


def generate_expressions(build_dir, functions, batches):
    import generate_expression_tests

    output = Path(build_dir).resolve() / "expression-sources"
    output.mkdir(parents=True, exist_ok=True)
    for source in output.glob("*test.cpp"):
        source.unlink()
    generate_expression_tests.src_folder = str(output) + "/"
    generate_expression_tests.main(functions, batches)
    if not list(output.glob("*test.cpp")):
        raise RuntimeError("Expression generator produced no tests")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("changed-unit-filter")
    distributions = commands.add_parser("changed-distribution-filter")
    distributions.add_argument("--all", action="store_true")
    expressions = commands.add_parser("generate-expressions")
    expressions.add_argument("build_dir")
    expressions.add_argument("--functions", nargs="*", default=[])
    expressions.add_argument("--batches", type=int, default=6)
    args = parser.parse_args()
    if args.command == "changed-unit-filter":
        print(changed_units())
    elif args.command == "changed-distribution-filter":
        print(changed_distributions([] if args.all else changed_files(), args.all))
    else:
        if args.batches < 1:
            parser.error("--batches must be positive")
        generate_expressions(args.build_dir, args.functions, args.batches)


if __name__ == "__main__":
    main()

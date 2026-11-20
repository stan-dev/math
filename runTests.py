#!/usr/bin/env python3
"""Build and run Stan Math tests with CMake."""

from concurrent.futures import ThreadPoolExecutor
import glob
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from argparse import ArgumentParser, RawTextHelpFormatter

def processCLIArgs(arguments):
    """
    Define and process the command line interface to the runTests.py script.
    """
    cli_description = "Generate and run stan math library tests."
    cli_epilog = "See more information at: https://github.com/stan-dev/math"

    parser = ArgumentParser(
        description=cli_description,
        epilog=cli_epilog,
        formatter_class=RawTextHelpFormatter,
    )

    parser.set_defaults(command=None)
    if arguments and arguments[0] in ("changed-unit-filter", "generate-expressions", "build-tests"):
        commands = parser.add_subparsers(dest="command", required=True)
        commands.add_parser("changed-unit-filter")
        expressions = commands.add_parser("generate-expressions")
        expressions.add_argument("build_dir")
        expressions.add_argument("--functions", nargs="*", default=[])
        expressions.add_argument("--batches", type=int, default=6)
        build = commands.add_parser("build-tests")
        build.add_argument("build_dir")
        build.add_argument("--source-regex", required=True)
        build.add_argument("--jobs", type=int, default=1)
        inputs = parser.parse_args(arguments)
        if getattr(inputs, "jobs", 1) < 1 or getattr(inputs, "batches", 1) < 1:
            parser.error("job and batch counts must be positive")
        return inputs

    parser.add_argument(
        "-j", metavar="N", type=int, default=1, help="parallel build jobs (default: 1)"
    )
    parser.add_argument(
        "--test-jobs", metavar="N", type=int, default=1,
        help="test executables run at once (default: 1)"
    )

    parser.add_argument(
        "-e",
        metavar="M",
        type=int,
        default=-1,
        help="number of files to split expressions tests in",
    )
    f_help_msg = "Only tests with file names matching these will be executed.\n"
    f_help_msg += "Example: '-f chol', '-f opencl', '-f prim'"
    parser.add_argument("-f", type=str, default=[], action="append", help=f_help_msg)
    changed_help = (
        "Use git to determine which tests may have changed, and run only those"
    )
    parser.add_argument("--changed", action="store_true", help=changed_help)
    parser.add_argument(
        "-d",
        "--debug",
        dest="debug",
        action="store_true",
        help="request additional script debugging output.",
    )
    parser.add_argument(
        "-m",
        "--make-only",
        dest="make_only",
        action="store_true",
        help="Don't run tests, just try to make them.",
    )
    parser.add_argument(
        "-t",
        "--test-only",
        dest="test_only",
        action="store_true",
        help="Don't make tests, just try to run them (must already be made).",
    )
    parser.add_argument(
        "--configure-only",
        dest="configure_only",
        action="store_true",
        help="Only configure the CMake build directory.",
    )
    parser.add_argument(
        "--run-all",
        dest="run_all",
        action="store_true",
        help="Don't stop at the first test failure, run all of them.",
    )
    parser.add_argument(
        "--only-functions",
        nargs="+",
        type=str,
        default=[],
        help="Function names to run expression tests for. Default: all functions",
    )
    parser.add_argument(
        "--jumbo-test",
        dest="do_jumbo",
        action="store_true",
        help="accepted for compatibility, directory tests already use one executable",
    )

    tests_help_msg = "The path(s) to the test case(s) to run.\n"
    tests_help_msg += "Example: 'test/unit', 'test/prob', and/or\n"
    tests_help_msg += "         'test/unit/math/prim/fun/abs_test.cpp', and/or\n"
    tests_help_msg += "         'test/prob/normal/normal_test.hpp' (a distribution fixture)"
    parser.add_argument("tests", nargs="*", type=str, help=tests_help_msg)
    parser.add_argument("--build-dir", default="build", help="CMake build directory")
    parser.add_argument("--fresh", action="store_true",
                        help="reset the CMake cache before configuring, like cmake --fresh")
    parser.add_argument("--cmake-arg", action="append", default=[],
                        help="extra configure argument, e.g. --cmake-arg=-DSTAN_THREADS=ON")
    parser.add_argument("--config", default="ReleaseTesting",
                        help="CMake build configuration (ReleaseTesting: -O1 without NDEBUG)")
    # And parse the command line against those rules
    inputs = parser.parse_args(arguments)
    if inputs.j < 1 or inputs.test_jobs < 1:
        parser.error("-j and --test-jobs must be positive")
    if inputs.test_only and (inputs.make_only or inputs.configure_only):
        parser.error("--test-only cannot be combined with --make-only or --configure-only")
    if inputs.fresh and inputs.test_only:
        parser.error("--fresh cannot be combined with --test-only")
    return inputs


def findChangedTests(files, existing_tests, limit=20):
    """Select priority tests from supplied changes and existing test paths."""
    dependencies = {"prim": ("prim", "rev", "mix", "fwd"),
                    "rev": ("rev", "mix"), "fwd": ("fwd", "mix")}
    selected = set()
    for path in files:
        if "opencl" in path:
            continue
        if path.endswith("_test.cpp"):
            selected.add(path)
        elif path.startswith("stan/") and path.endswith(".hpp"):
            candidate = path.replace("stan/", "test/unit/", 1)[:-4] + "_test.cpp"
            layer = candidate.split("/")[3]
            for replacement in dependencies.get(layer, (layer,)):
                test = candidate.replace("/" + layer + "/", "/" + replacement + "/", 1)
                if test in existing_tests:
                    selected.add(test)
    return sorted(selected) if len(selected) <= limit else []


def is_fixture(path):
    return re.fullmatch(r"test/prob/[^/]+/[^/]+_test\.hpp", path) is not None


def cmake_escape(text):
    return re.sub(r"([.+*?^$()\[\]{}|\\])", r"\\\1", text)


def selection_regex(paths):
    """Match the source paths of the selected files and fixtures' generated tests."""
    patterns = [cmake_escape(path[:-len("_test.hpp")]) + "_[0-9]+_generated_.*" if is_fixture(path)
                else cmake_escape(path) for path in sorted(paths)
                if is_fixture(path) or path.endswith(".cpp")]
    return "^(" + "|".join(patterns) + ")$" if patterns else ""


def fixture_tests(inventory, fixtures):
    """Expand fixtures to the generated sources CMake registered for them."""
    prefixes = tuple(fixture[:-len("_test.hpp")] + "_" for fixture in fixtures)
    sources = set()
    for test in inventory["tests"]:
        for entry in test.get("properties", []):
            if entry["name"] != "LABELS":
                continue
            for label in entry["value"]:
                for prefix in prefixes:
                    if label.startswith(prefix) and re.match(
                            r"[0-9]+_generated_.*_test\.cpp$", label[len(prefix):]):
                        sources.add(label)
    return sorted(sources)


def select_paths(matches, directories):
    """Validate expanded paths and remove selections covered by a directory."""
    paths = set()
    for argument, expanded in matches.items():
        if not expanded:
            raise ValueError("No matching tests: " + argument)
        for path in expanded:
            if not any(path == base or path.startswith(base + "/")
                       for base in ("test/unit", "test/prob", "test/expressions")):
                raise ValueError("Unsupported test path: " + path)
            if path not in directories and not path.endswith("_test.cpp") and not is_fixture(path):
                raise ValueError("Expected a test directory, *_test.cpp, or "
                                 "test/prob/*/*_test.hpp: " + path)
            paths.add(path)
    return [path for path in sorted(paths)
            if not any(path.startswith(parent + "/") for parent in directories)]


def target_for(path):
    if path.endswith(".cpp"):
        target = path[:-4].replace("/", "_")
    else:
        prefix = path.replace("test/unit/math", "unit_math", 1)
        prefix = prefix.replace("test/unit", "unit", 1)
        prefix = prefix.replace("test/prob", "prob_math", 1)
        prefix = prefix.replace("test/expressions", "expressions", 1)
        target = prefix.replace("/", "_") + "_all_tests"
    folder = "test/unit" if path.startswith("test/unit") else (
        "test/prob" if path.startswith("test/prob") else "test")
    return target, folder


def build_commands(inputs, paths, root, build_dir):
    """Plan configuration and builds for the CMake targets selected by the user."""
    if inputs.test_only:
        return []
    expressions = any(path.startswith("test/expressions") for path in paths)
    regex = selection_regex(paths)
    configure = ["cmake", "-S", str(root), "-B", str(build_dir),
                 "-DCMAKE_BUILD_TYPE=" + inputs.config, "-DBUILD_TESTING=ON",
                 "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON", "-DSTAN_TEST_HEADERS=OFF",
                 "-DSTAN_BUILD_PROB_TESTS=" + ("ON" if any(path.startswith("test/prob") for path in paths) else "OFF"),
                 "-DSTAN_BUILD_EXPRESSION_TESTS=" + ("ON" if expressions else "OFF"),
                 "-DSTAN_TEST_NAME_FILTERS=" + ";".join(inputs.f),
                 "-DSTAN_SELECTED_TEST_REGEX=" + regex]
    if expressions:
        configure.append("-DSTAN_EXPRESSION_TEST_DIR=" + str(build_dir / "expression-sources"))
    if inputs.configure_only:
        return [configure + inputs.cmake_arg]
    # Individual files and fixtures build through one aggregate target.
    targets = [target_for(path)[0] for path in paths if not path.endswith(".cpp") and not is_fixture(path)]
    if regex:
        targets.append("stan_selected_tests")
    build = ["cmake", "--build", str(build_dir), "--config", inputs.config,
             "--parallel", str(inputs.j), "--target"] + targets
    return [configure + inputs.cmake_arg, build]


def parse_cache(text):
    cache = {}
    for line in text.splitlines():
        if not line.startswith(("//", "#")) and "=" in line and ":" in line.split("=", 1)[0]:
            key, value = line.split("=", 1)
            cache[key.split(":", 1)[0]] = value
    return cache


def test_command(path, build_dir, config, cache, windows):
    """Plan an executable invocation using supplied CMake settings."""
    target, folder = target_for(path)
    executable = build_dir / folder / target
    if cache.get("CMAKE_CONFIGURATION_TYPES"):
        executable = executable.parent / config / executable.name
    if cache.get("CMAKE_SYSTEM_NAME") == "Windows" or windows:
        executable = executable.with_suffix(".exe")
    args = [str(executable), "--gtest_output=xml:" + str(build_dir / "test-results" / (target + ".xml"))]
    if cache.get("STAN_EFFECTIVE_MPI", cache.get("STAN_MPI")) in ("ON", "TRUE", "1"):
        ranks = "2" if not path.endswith(".cpp") or "mpi_" in path else "1"
        args = [cache.get("MPIEXEC_EXECUTABLE", "mpirun"),
                cache.get("MPIEXEC_NUMPROC_FLAG", "-np"), ranks] + [
                    flag for flag in cache.get("MPIEXEC_PREFLAGS", "").split(";") if flag] + args + [
                    flag for flag in cache.get("MPIEXEC_POSTFLAGS", "").split(";") if flag]
    return args


def source_filter(paths):
    return "^(" + "|".join(re.escape(path) for path in sorted(set(paths))) + ")$" if paths else ""


def build_tests(inventory, source_pattern, build_dir, jobs, batch_size):
    """Plan builds using the source labels registered by CMake."""
    pattern = re.compile(source_pattern)
    selected, groups = [], {}
    for test in inventory["tests"]:
        labels = [label for entry in test.get("properties", [])
                  if entry["name"] == "LABELS" for label in entry["value"]]
        if any(pattern.search(label) for label in labels if not label.startswith("stan_build:")):
            selected.append(test["name"])
        for label in labels:
            if label.startswith("stan_build:"):
                groups.setdefault(label[len("stan_build:"):], set()).add(test["name"])
    if not selected:
        raise ValueError("No registered tests match: " + source_pattern)
    # Makefile generators serialize separate top-level targets. An existing
    # aggregate shares one dependency graph and honors the requested job count.
    remaining, targets = set(selected), []
    for target, members in sorted(groups.items(), key=lambda item: (-len(item[1]), item[0])):
        if len(members) > 1 and members <= remaining:
            targets.append(target)
            remaining.difference_update(members)
    targets.extend(name for name in selected if name in remaining)
    return [["cmake", "--build", str(build_dir), "--parallel", str(jobs), "--target"]
            + targets[offset:offset + batch_size] for offset in range(0, len(targets), batch_size)]


def main(arguments, root, working_dir, windows, error_stream):
    """Perform command-line I/O. Helpers above only transform explicit inputs."""
    inputs = processCLIArgs(arguments)
    try:
        if getattr(inputs, "changed", False) or inputs.command == "changed-unit-filter":
            files = subprocess.run(["git", "diff", "--name-only", "--diff-filter=d",
                                    "origin/develop...HEAD"], cwd=root, check=True,
                                   text=True, capture_output=True).stdout.splitlines()
            if getattr(inputs, "debug", False):
                print("Changed files:", files)
            existing = {path.relative_to(root).as_posix()
                        for path in (root / "test/unit").rglob("*_test.cpp")}
            paths = findChangedTests(files, existing)
            if inputs.command == "changed-unit-filter":
                print(source_filter(paths))
                return 0
            if not paths:
                print("No changed tests were found!")
                return 0
        elif inputs.command is None:
            matches = {}
            for argument in inputs.tests:
                # os.path.join keeps a trailing slash, so "test/prob/*/" matches directories.
                expanded = glob.glob(os.path.join(str(root), argument))
                if not expanded and argument.replace("\\", "/").lstrip("./").startswith(
                        ("test/prob/", "test/expressions/")) and argument.endswith("_test.cpp"):
                    expanded = [str(root / argument)]
                matches[argument] = [Path(path).resolve().relative_to(root).as_posix()
                                     for path in expanded]
            directories = {path for expanded in matches.values() for path in expanded
                           if (root / path).is_dir()}
            paths = select_paths(matches, directories)
        build_dir = (working_dir / inputs.build_dir).resolve()
        if inputs.command == "build-tests":
            inventory = subprocess.run(["ctest", "--show-only=json-v1"],
                                       cwd=build_dir, check=True, text=True, capture_output=True)
            commands = build_tests(json.loads(inventory.stdout), inputs.source_regex,
                                   build_dir, inputs.jobs, 50 if windows else 500)
        else:
            generate = inputs.command == "generate-expressions"
            if not generate:
                paths = [path for path in paths if not path.endswith(".cpp")
                         or all(value in path for value in inputs.f)]
                if not paths and not inputs.configure_only:
                    raise ValueError("No matching tests found")
                expressions = any(path.startswith("test/expressions") for path in paths)
                if inputs.only_functions and not expressions:
                    raise ValueError("--only-functions requires test/expressions")
                generate = expressions and not inputs.test_only
            if generate:
                # Keep the legacy generator's output setting inside the I/O boundary.
                sys.path.insert(0, str(root / "test"))
                import generate_expression_tests
                output = build_dir / "expression-sources"
                output.mkdir(parents=True, exist_ok=True)
                for source in output.glob("*test.cpp"):
                    source.unlink()
                functions = inputs.functions if inputs.command else inputs.only_functions
                # Like develop, four files per job keep each translation unit's memory bounded.
                batches = inputs.batches if inputs.command else (inputs.e if inputs.e > 0 else inputs.j * 4)
                generate_expression_tests.src_folder = str(output) + "/"
                generate_expression_tests.main(functions, batches)
                if not list(output.glob("*test.cpp")):
                    raise RuntimeError("Expression generator produced no tests")
            if inputs.command == "generate-expressions":
                return 0
            commands = build_commands(inputs, paths, root, build_dir)
            if inputs.fresh:
                # cmake --fresh (CMake 3.24+) deletes the same two paths.
                (build_dir / "CMakeCache.txt").unlink(missing_ok=True)
                shutil.rmtree(build_dir / "CMakeFiles", ignore_errors=True)
        for command in commands:
            print(" ".join(shlex.quote(str(arg)) for arg in command), flush=True)
            subprocess.run(command, cwd=root, check=True)
        if inputs.command == "build-tests" or inputs.make_only or inputs.configure_only:
            return 0
        fixtures = [path for path in paths if is_fixture(path)]
        if fixtures:
            inventory = subprocess.run(["ctest", "--show-only=json-v1"],
                                       cwd=build_dir, check=True, text=True, capture_output=True)
            generated =fixture_tests(json.loads(inventory.stdout), fixtures)
            if not generated:
                raise ValueError("No generated tests found for: " + " ".join(fixtures))
            paths = [path for path in paths if not is_fixture(path)] + generated
        cache_file = build_dir / "CMakeCache.txt"
        cache = parse_cache(cache_file.read_text() if cache_file.exists() else "")
        (build_dir / "test-results").mkdir(parents=True, exist_ok=True)
        failed = False
        commands = [test_command(path, build_dir, inputs.config, cache, windows) for path in paths]
        if inputs.test_jobs > 1:
            # Concurrent executables print their output in order once each finishes.
            def run(command):
                return subprocess.run(command, cwd=root, text=True,
                                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            with ThreadPoolExecutor(inputs.test_jobs) as pool:
                for command, result in zip(commands, pool.map(run, commands)):
                    print(" ".join(shlex.quote(arg) for arg in command), flush=True)
                    print(result.stdout, end="", flush=True)
                    failed = failed or result.returncode != 0
            return int(failed)
        for command in commands:
            print(" ".join(shlex.quote(arg) for arg in command), flush=True)
            try:
                subprocess.run(command, cwd=root, check=True)
            except subprocess.CalledProcessError:
                if not inputs.run_all:
                    raise
                failed = True
        return int(failed)
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(error, file=error_stream)
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:], Path(__file__).resolve().parent,
                  Path.cwd(), os.name == "nt", sys.stderr))

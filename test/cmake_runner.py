"""CMake backend for runTests.py; directories select one grouped executable."""

import glob
import os
from pathlib import Path
import re
import shlex
import subprocess

ROOT = Path(__file__).resolve().parent.parent


def select_paths(arguments):
    paths = set()
    for argument in arguments:
        matches = glob.glob(str(ROOT / argument))
        # Distribution sources are generated in the build tree, not the source tree.
        if not matches and argument.replace("\\", "/").lstrip("./").startswith(("test/prob/", "test/expressions/")) and argument.endswith("_test.cpp"):
            matches = [str(ROOT / argument)]
        if not matches:
            raise ValueError("No matching tests: " + argument)
        for match in matches:
            path = Path(match).resolve().relative_to(ROOT).as_posix()
            if not (path == "test/unit" or path.startswith("test/unit/")
                    or path == "test/prob" or path.startswith("test/prob/")
                    or path == "test/expressions" or path.startswith("test/expressions/")):
                raise ValueError("Unsupported test path: " + path)
            if not (ROOT / path).is_dir() and not path.endswith("_test.cpp"):
                raise ValueError("Expected a test directory or *_test.cpp: " + path)
            paths.add(path)
    # A parent directory already covers its descendants.
    return [path for path in sorted(paths)
            if not any(path.startswith(parent + "/") for parent in paths
                       if (ROOT / parent).is_dir())]


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


def selection_filter(paths, filters):
    parts = []
    for path in paths:
        if path.endswith(".cpp"):
            if all(value in path for value in filters):
                parts.append(re.escape(path) + "$")
        else:
            parts.append(re.escape(path) + "/")
    if not parts:
        raise ValueError("No matching tests found")
    return "^(" + "|".join(parts) + ")"


def command(args):
    print(" ".join(shlex.quote(str(arg)) for arg in args), flush=True)
    subprocess.run(args, cwd=ROOT, check=True)


def read_cache(build_dir):
    cache = {}
    cache_file = build_dir / "CMakeCache.txt"
    if cache_file.exists():
        for line in cache_file.read_text().splitlines():
            if not line.startswith(("//", "#")) and "=" in line and ":" in line.split("=", 1)[0]:
                key, value = line.split("=", 1)
                cache[key.split(":", 1)[0]] = value
    return cache


def run(inputs):
    arguments = inputs.tests
    if inputs.changed:
        from runTests import findChangedTests
        arguments = [re.sub(r"\.exe$", "", path) + ".cpp" for path in findChangedTests(inputs.debug)]
        if not arguments:
            print("No changed tests were found!")
            return 0
    paths = select_paths(arguments)
    if not paths:
        raise ValueError("Specify a test directory or source file")
    if inputs.f:
        paths = [path for path in paths if not path.endswith(".cpp")
                 or all(value in path for value in inputs.f)]
        if not paths:
            raise ValueError("No matching tests found")
    expressions = any(path.startswith("test/expressions") for path in paths)
    if inputs.only_functions and not expressions:
        raise ValueError("--only-functions requires test/expressions")
    build_dir = Path(inputs.build_dir).resolve()
    if not inputs.test_only:
        if expressions:
            from cmake_ci import generate_expressions
            generate_expressions(build_dir, inputs.only_functions,
                                 inputs.e if inputs.e > 0 else inputs.j)
        configure = ["cmake", "-S", str(ROOT), "-B", str(build_dir),
                     "-DCMAKE_BUILD_TYPE=" + inputs.config, "-DBUILD_TESTING=ON",
                     "-DSTAN_TEST_HEADERS=OFF",
                     "-DSTAN_BUILD_PROB_TESTS=" + ("ON" if any(path.startswith("test/prob") for path in paths) else "OFF"),
                     "-DSTAN_BUILD_EXPRESSION_TESTS=" + ("ON" if expressions else "OFF"),
                     # Keep ordinary directory/file targets configured so switching
                     # selections preserves their generated build metadata.
                     "-DSTAN_TEST_FILTER=" + (selection_filter(paths, inputs.f) if inputs.f else ""),
                     "-DSTAN_TEST_NAME_FILTERS=" + ";".join(inputs.f)]
        if expressions:
            configure.append("-DSTAN_EXPRESSION_TEST_DIR=" + str(build_dir / "expression-sources"))
        command(configure + inputs.cmake_arg)
        command(["cmake", "--build", str(build_dir), "--config", inputs.config,
                 "--parallel", str(inputs.j), "--target"] + [target_for(path)[0] for path in paths])
    if inputs.make_only:
        return 0
    cache = read_cache(build_dir)
    results = build_dir / "test-results"
    results.mkdir(parents=True, exist_ok=True)
    failed = False
    for path in paths:
        target, folder = target_for(path)
        executable = build_dir / folder / target
        if cache.get("CMAKE_CONFIGURATION_TYPES"):
            executable = executable.parent / inputs.config / executable.name
        if cache.get("CMAKE_SYSTEM_NAME") == "Windows" or os.name == "nt":
            executable = executable.with_suffix(".exe")
        args = [str(executable), "--gtest_output=xml:" + str(results / (target + ".xml"))]
        if cache.get("STAN_MPI") in ("ON", "TRUE", "1"):
            # Grouped MPI suites can contain collective tests requiring two ranks.
            ranks = "2" if not path.endswith(".cpp") or "mpi_" in path else "1"
            args = [cache.get("MPIEXEC_EXECUTABLE", "mpirun"),
                    cache.get("MPIEXEC_NUMPROC_FLAG", "-np"), ranks] + [
                        flag for flag in cache.get("MPIEXEC_PREFLAGS", "").split(";") if flag] + args + [
                        flag for flag in cache.get("MPIEXEC_POSTFLAGS", "").split(";") if flag]
        try:
            command(args)
        except subprocess.CalledProcessError:
            if not inputs.run_all:
                raise
            failed = True
    return int(failed)

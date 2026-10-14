#!/usr/bin/env python3
"""CMake test runner, CI commands, and regression checks.

Run regression checks with: python3 test/cmake.py self-test
"""

import argparse
import glob
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent.parent
MODULE = sys.modules[__name__]

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
        arguments = findChangedTests(inputs.debug)
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
            generate_expressions(build_dir, inputs.only_functions,
                                 inputs.e if inputs.e > 0 else inputs.j)
        configure = ["cmake", "-S", str(ROOT), "-B", str(build_dir),
                     "-DCMAKE_BUILD_TYPE=" + inputs.config, "-DBUILD_TESTING=ON",
                     "-DSTAN_TEST_HEADERS=OFF",
                     "-DSTAN_BUILD_PROB_TESTS=" + ("ON" if any(path.startswith("test/prob") for path in paths) else "OFF"),
                     "-DSTAN_BUILD_EXPRESSION_TESTS=" + ("ON" if expressions else "OFF"),
                     "-DSTAN_TEST_NAME_FILTERS=" + ";".join(inputs.f)]
        if expressions:
            configure.append("-DSTAN_EXPRESSION_TEST_DIR=" + str(build_dir / "expression-sources"))
        if inputs.fresh:
            configure.append("--fresh")
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
        if cache.get("STAN_EFFECTIVE_MPI", cache.get("STAN_MPI")) in ("ON", "TRUE", "1"):
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

    return source_filter(findChangedTests(False))


def build_tests(build_dir, source_pattern, jobs):
    inventory = subprocess.run(["ctest", "--test-dir", str(build_dir), "--show-only=json-v1"],
                               check=True, text=True, capture_output=True)
    selected = []
    pattern = re.compile(source_pattern)
    for test in json.loads(inventory.stdout)["tests"]:
        labels = next((entry["value"] for entry in test.get("properties", [])
                       if entry["name"] == "LABELS"), [])
        if any(pattern.search(label) for label in labels):
            selected.append(test["name"])
    if not selected:
        raise ValueError("No registered tests match: " + source_pattern)
    batch_size = 50 if sys.platform == "win32" else 500
    for offset in range(0, len(selected), batch_size):
        subprocess.run(["cmake", "--build", str(build_dir), "--parallel", str(jobs),
                        "--target"] + selected[offset:offset + batch_size], check=True)


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
    commands.add_parser("self-test")
    commands.add_parser("changed-unit-filter")
    distributions = commands.add_parser("changed-distribution-filter")
    distributions.add_argument("--all", action="store_true")
    expressions = commands.add_parser("generate-expressions")
    expressions.add_argument("build_dir")
    expressions.add_argument("--functions", nargs="*", default=[])
    expressions.add_argument("--batches", type=int, default=6)
    build = commands.add_parser("build-tests")
    build.add_argument("build_dir")
    build.add_argument("--source-regex", required=True)
    build.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    if args.command == "self-test":
        global runTests
        sys.path.insert(0, str(ROOT))
        import runTests
        unittest.main(argv=[sys.argv[0]])
    elif args.command == "changed-unit-filter":
        print(changed_units())
    elif args.command == "changed-distribution-filter":
        print(changed_distributions([] if args.all else changed_files(), args.all))
    elif args.command == "build-tests":
        if args.jobs < 1:
            parser.error("--jobs must be positive")
        try:
            build_tests(args.build_dir, args.source_regex, args.jobs)
        except ValueError as error:
            parser.error(str(error))
    else:
        if args.batches < 1:
            parser.error("--batches must be positive")
        generate_expressions(args.build_dir, args.functions, args.batches)


class CMakeRunnerTest(unittest.TestCase):
    def inputs(self, build_dir, **overrides):
        values = dict(tests=["./test/unit/math/prim/fun"], changed=False,
                      debug=False, f=[], only_functions=[], build_dir=build_dir,
                      test_only=False, make_only=False, config="Release", e=-1,
                      j=6, cmake_arg=[], run_all=False, fresh=False)
        values.update(overrides)
        return SimpleNamespace(**values)

    def test_directory_builds_and_runs_one_group(self):
        with tempfile.TemporaryDirectory() as build_dir, patch.object(MODULE, "command") as command:
            self.assertEqual(run(self.inputs(build_dir)), 0)
            configure, build, execute = [call.args[0] for call in command.call_args_list]
            self.assertEqual(build[build.index("--target") + 1:], ["unit_math_prim_fun_all_tests"])
            self.assertEqual(Path(execute[0]).name, "unit_math_prim_fun_all_tests")
            self.assertIn("test-results/unit_math_prim_fun_all_tests.xml", execute[1])
            self.assertNotIn("--fresh", configure)


    @unittest.skipUnless(shutil.which("cmake"), "CMake is required")
    def test_fresh_resets_cached_settings_and_preserves_build_files(self):
        version = subprocess.run(["cmake", "--version"], check=True,
                                 text=True, stdout=subprocess.PIPE).stdout
        match = re.search(r"cmake version (\d+)\.(\d+)", version)
        if tuple(map(int, match.groups())) < (3, 24):
            self.skipTest("CMake --fresh requires CMake 3.24")
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source"
            source.mkdir()
            fixture = source / "test/unit/math/prim/fun/sum_test.cpp"
            fixture.parent.mkdir(parents=True)
            fixture.touch()
            (source / "CMakeLists.txt").write_text("""
cmake_minimum_required(VERSION 3.24)
project(fresh_fixture NONE)
option(STAN_THREADS "Thread support" OFF)
add_custom_target(test_unit_math_prim_fun_sum_test
  COMMAND "${CMAKE_COMMAND}" -E touch "${CMAKE_BINARY_DIR}/built-marker")
""")
            build = Path(temporary) / "build"
            subprocess.run(["cmake", "-S", str(source), "-B", str(build),
                            "-DSTAN_THREADS=ON"], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            (build / "keep-marker").touch()
            (build / "CMakeFiles/remove-marker").touch()
            with patch.object(MODULE, "ROOT", source):
                run(self.inputs(str(build), make_only=True,
                    tests=["test/unit/math/prim/fun/sum_test.cpp"]))
                self.assertEqual(read_cache(build)["STAN_THREADS"], "ON")
                with patch.object(sys, "argv", ["runTests.py", "--fresh", "--make-only",
                        "--build-dir", str(build), "test/unit/math/prim/fun/sum_test.cpp"]):
                    run(runTests.processCLIArgs())
            cache = read_cache(build)
            self.assertEqual(cache["STAN_THREADS"], "OFF")
            self.assertTrue((build / "keep-marker").exists())
            self.assertTrue((build / "built-marker").exists())
            self.assertFalse((build / "CMakeFiles/remove-marker").exists())


    def test_overlap_runs_only_parent(self):
        self.assertEqual(select_paths([
            "test/unit/math/prim/fun/sum_test.cpp", "./test/unit/math/prim/fun",
            "test/unit/math/prim", "test/unit/math/prim"]), ["test/unit/math/prim"])

    def test_filters_do_not_turn_directory_into_individual_targets(self):
        with tempfile.TemporaryDirectory() as build_dir, patch.object(MODULE, "command") as command:
            run(self.inputs(build_dir, f=["sum", "prim"], make_only=True))
            self.assertIn("-DSTAN_TEST_NAME_FILTERS=sum;prim", command.call_args_list[0].args[0])
            self.assertEqual(command.call_args_list[1].args[0][-1], "unit_math_prim_fun_all_tests")
            self.assertEqual(command.call_count, 2)

    def test_test_only_skips_build(self):
        with tempfile.TemporaryDirectory() as build_dir, patch.object(MODULE, "command") as command:
            run(self.inputs(build_dir, test_only=True))
            self.assertEqual(command.call_count, 1)
            self.assertEqual(Path(command.call_args.args[0][0]).name, "unit_math_prim_fun_all_tests")

    def test_run_all_preserves_failure_and_continues(self):
        with tempfile.TemporaryDirectory() as build_dir, patch.object(MODULE, "command") as command:
            command.side_effect = [subprocess.CalledProcessError(1, "failed-test"), None]
            result = run(self.inputs(build_dir, tests=[
                "test/unit/math/prim/fun/sum_test.cpp",
                "test/unit/math/prim/fun/abs_test.cpp"], test_only=True, run_all=True))
            self.assertEqual(result, 1)
            self.assertEqual(command.call_count, 2)

    def test_build_failure_stops_before_execution(self):
        with tempfile.TemporaryDirectory() as build_dir, patch.object(MODULE, "command") as command:
            command.side_effect = [None, subprocess.CalledProcessError(2, "build")]
            with self.assertRaises(subprocess.CalledProcessError):
                run(self.inputs(build_dir))
            self.assertEqual(command.call_count, 2)


    def test_cli_jobs_reach_build_command(self):
        for flags, jobs in (([], 1), (["-j20"], 20)):
            with self.subTest(jobs=jobs), tempfile.TemporaryDirectory() as directory, patch.object(
                    sys, "argv", ["runTests.py", "--make-only", "--build-dir", directory] + flags + [
                    "test/unit/math/prim/fun"]), patch.object(MODULE, "command") as command:
                inputs = runTests.processCLIArgs()
                self.assertFalse(inputs.fresh)
                run(inputs)
                build = command.call_args.args[0]
                self.assertEqual(build[build.index("--parallel") + 1], str(jobs))

    def test_cli_rejects_invalid_options(self):
        for flags in (["-j0"], ["--fresh", "--test-only"], ["--backend", "make"]):
            with self.subTest(flags=flags), patch.object(sys, "argv", ["runTests.py"] + flags + [
                    "test/unit"]), patch.object(sys, "stderr"):
                with self.assertRaises(SystemExit) as error:
                    runTests.processCLIArgs()
                self.assertEqual(error.exception.code, 2)

    def test_invalid_path_rejected(self):
        with self.assertRaises(ValueError):
            select_paths(["stan/math"])

    def test_filtered_file_selection_cannot_succeed_without_tests(self):
        with self.assertRaises(ValueError):
            run(self.inputs("build", tests=[
                "test/unit/math/prim/fun/sum_test.cpp"], f=["cholesky"], test_only=True))

    def test_multi_configuration_executable_path(self):
        with tempfile.TemporaryDirectory() as build_dir, patch.object(MODULE, "command") as command:
            Path(build_dir, "CMakeCache.txt").write_text("CMAKE_CONFIGURATION_TYPES:STRING=Debug;Release\n")
            run(self.inputs(build_dir, test_only=True, config="Debug"))
            self.assertEqual(Path(command.call_args.args[0][0]),
                             Path(build_dir, "test/unit/Debug/unit_math_prim_fun_all_tests"))

    def test_effective_mpi_setting_controls_launch(self):
        with tempfile.TemporaryDirectory() as build_dir, patch.object(MODULE, "command") as command:
            Path(build_dir, "CMakeCache.txt").write_text(
                "STAN_EFFECTIVE_MPI:INTERNAL=ON\nMPIEXEC_EXECUTABLE:FILEPATH=mpirun\n")
            run(self.inputs(build_dir, test_only=True))
            self.assertEqual(command.call_args.args[0][:3], ["mpirun", "-np", "2"])


class CmakeCISelectionTest(unittest.TestCase):
    def test_unit_only_changes_skip_distributions(self):
        self.assertEqual(changed_distributions(
            ["test/unit/math/prim/fun/sum_test.cpp"]), "")

    def test_fixture_change_selects_all_modes_and_batches(self):
        selected = changed_distributions(
            ["test/prob/bernoulli/bernoulli_test.hpp"])
        for mode in ("v", "ffv", "vv", "fd", "fv", "ffd"):
            for batch in ("00000", "00010"):
                self.assertRegex(
                    f"test/prob/bernoulli/bernoulli_{batch}_generated_{mode}_test.cpp", selected)
        self.assertIsNone(re.search(selected,
            "test/prob/bernoulli/bernoulli_cdf_00000_generated_v_test.cpp"))

    def test_shared_dependencies_select_all_distributions(self):
        for dependency in ("stan/math/prim/fun/sum.hpp", "lib/boost_1.87.0/boost/config.hpp",
                           "test/prob/utility.hpp", "test/unit/math/expect_near_rel.hpp",
                           "test/prob/generate_tests.cpp", "cmake/CMakeLists.txt", "Jenkinsfile"):
            with self.subTest(dependency=dependency):
                selected = changed_distributions([dependency])
                self.assertRegex("test/prob/bernoulli/bernoulli_00000_generated_v_test.cpp", selected)
                self.assertRegex("test/prob/normal/normal_00000_generated_ffv_test.cpp", selected)


    def test_build_selection_uses_registered_labels(self):
        inventory = {"tests": [
            {"name": "test_unit_math_prim_fun_sum_test", "properties": [
                {"name": "LABELS", "value": ["prim", "test/unit/math/prim/fun/sum_test.cpp"]}]},
            {"name": "test_unit_math_rev_core_grad_test", "properties": [
                {"name": "LABELS", "value": ["rev", "test/unit/math/rev/core/grad_test.cpp"]}]},
        ]}
        with patch.object(subprocess, "run") as command:
            command.return_value = subprocess.CompletedProcess([], 0, json.dumps(inventory))
            build_tests("build", "^test/unit/math/prim/", 12)
            self.assertEqual(command.call_args.args[0], [
                "cmake", "--build", "build", "--parallel", "12", "--target",
                "test_unit_math_prim_fun_sum_test"])

    def test_empty_build_selection_fails(self):
        with patch.object(subprocess, "run") as command:
            command.return_value = subprocess.CompletedProcess([], 0, '{"tests": []}')
            with self.assertRaisesRegex(ValueError, "No registered tests match"):
                build_tests("build", "^test/unit/", 1)


class CMakeBuildTest(unittest.TestCase):
    def check_build(self, pch_enabled):
        with tempfile.TemporaryDirectory(prefix="stan-cmake-build-") as temporary:
            source = Path(temporary) / "source"
            build = Path(temporary) / "build"
            commands = Path(temporary) / "compile.jsonl"

            def write(path, text):
                destination = source / path
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_text(text)

            def run(*args):
                result = subprocess.run(args, text=True, stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT)
                self.assertEqual(result.returncode, 0, result.stdout)

            write("launcher.py", """import json, subprocess, sys
with open(sys.argv[1], 'a') as log:
    log.write(json.dumps(sys.argv[2:]) + '\\n')
sys.exit(subprocess.call(sys.argv[2:]))
""")
            write("CMakeLists.txt", f"""
cmake_minimum_required(VERSION 3.16.3)
project(stan_cmake_fixture LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(STAN_TEST_PCH {'ON' if pch_enabled else 'OFF'})
set(STAN_TEST_OPTIMIZATION_FLAGS -O1)
add_library(stanmath INTERFACE)
target_include_directories(stanmath INTERFACE "${{CMAKE_CURRENT_SOURCE_DIR}}")
add_library(gtest_main STATIC main.cpp)
enable_testing()
include("{(ROOT / 'test/CMakeLists.txt').as_posix()}")
""")
            write("cmake/pch.cpp", "// Empty executable/PCH source.\n")
            write("gtest/gtest.h", """#pragma once

#include <vector>

using test_callback = bool (*)();
inline std::vector<test_callback>& test_callbacks() {
  static std::vector<test_callback> callbacks;
  return callbacks;
}
struct test_registration {
  explicit test_registration(test_callback callback) {
    test_callbacks().push_back(callback);
  }
};
""")
            write("main.cpp", """#include <gtest/gtest.h>
int main() {
  for (auto callback : test_callbacks()) {
    if (!callback()) return 1;
  }
  return 0;
}
""")
            write("stan/math/prim.hpp", """#pragma once

#ifdef STAN_PCH_TEST_SETTING
constexpr int configured_value = STAN_PCH_TEST_SETTING;
#else
constexpr int configured_value = 1;
#endif
""")
            for layer in ("fwd", "rev", "mix"):
                write(f"stan/math/{layer}.hpp", "#include <stan/math/prim.hpp>\n")
            write("stan/math/mix.hpp", "#include <stan/math/prim.hpp>\nusing mix_only_alias = double;\n")
            write("test/unit/math/test_ad.hpp", """#pragma once

#include <stan/math/mix.hpp>
#include <gtest/gtest.h>

using test_alias = double;
""")
            write("unit/math/prim/fun/plain_test.cpp", """#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

using mix_only_alias = int;
using test_alias = int;
bool plain_test() { return configured_value == 1; }
test_registration plain_registration(plain_test);
""")
            write("unit/math/prim/fun/ad_test.cpp", """#include <test/unit/math/test_ad.hpp>
bool ad_test() { return configured_value == 1 && sizeof(test_alias) == sizeof(double); }
test_registration ad_registration(ad_test);
""")
            write("unit/math/prim/fun/configuration_test.cpp", """#define STAN_PCH_TEST_SETTING 2

#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

static_assert(configured_value == 2, "Configuration must precede Stan headers");
bool configuration_test() { return configured_value == 2; }
test_registration configuration_registration(configuration_test);
""")
            write("generated_template.cpp", """#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

bool generated_test() { return configured_value == 1; }
test_registration generated_registration(generated_test);
""")
            write("unit/CMakeLists.txt", """
foreach(layer prim fwd rev mix)
  stan_add_pch(unit_${layer}_pch <gtest/gtest.h> <stan/math/${layer}.hpp>)
endforeach()
stan_add_pch(unit_ad_pch <test/unit/math/test_ad.hpp>)
add_gtest_grouped_test("${CMAKE_CURRENT_SOURCE_DIR}/math" unit_math unit_mix_pch)
set(generated "${PROJECT_BINARY_DIR}/generated/generated_test.cpp")
add_custom_command(OUTPUT "${generated}"
  COMMAND "${CMAKE_COMMAND}" -E make_directory "${PROJECT_BINARY_DIR}/generated"
  COMMAND "${CMAKE_COMMAND}" -E copy "${PROJECT_SOURCE_DIR}/generated_template.cpp" "${generated}"
  DEPENDS "${PROJECT_SOURCE_DIR}/generated_template.cpp")
add_custom_target(generate_fixture DEPENDS "${generated}")
set_source_files_properties("${generated}" PROPERTIES STAN_GENERATOR_TARGET generate_fixture)
add_gtest_grouped_test("${PROJECT_BINARY_DIR}/generated" expressions unit_prim_pch "${generated}")
""")
            launcher = f"{sys.executable};{source / 'launcher.py'};{commands}"
            run("cmake", "-S", str(source), "-B", str(build),
                "-DCMAKE_CXX_COMPILER_LAUNCHER=" + launcher)
            run("cmake", "--build", str(build), "--parallel", "1", "--target",
                "unit_math_prim_fun_all_tests", "expressions_all_tests")
            run(str(build / "unit/unit_math_prim_fun_all_tests"))
            run(str(build / "unit/expressions_all_tests"))

            compiled = [json.loads(line) for line in commands.read_text().splitlines()]
            fixtures = {name: source / f"unit/math/prim/fun/{name}_test.cpp"
                        for name in ("plain", "ad", "configuration")}
            fixtures["generated"] = build / "generated/generated_test.cpp"
            expected_pch = {"plain": "unit_prim_pch", "ad": "unit_ad_pch",
                            "configuration": None, "generated": "unit_prim_pch"}
            for name, path in fixtures.items():
                invocations = [args for args in compiled if str(path) in args]
                self.assertEqual(len(invocations), 1, (name, invocations))
                pch = expected_pch[name] if pch_enabled else None
                pch_arguments = [arg for arg in invocations[0] if "cmake_pch" in arg]
                if pch:
                    self.assertTrue(any(pch in arg for arg in pch_arguments), invocations[0])
                else:
                    self.assertFalse(pch_arguments, invocations[0])

            # Switching executables must not compile any test source again.
            run("cmake", "--build", str(build), "--parallel", "1", "--target",
                "test_unit_math_prim_fun_plain_test", "unit_math_all_tests",
                "test_expressions_generated_test")
            run(str(build / "unit/test_unit_math_prim_fun_plain_test"))
            run(str(build / "unit/unit_math_all_tests"))
            run(str(build / "unit/test_expressions_generated_test"))
            after = [json.loads(line) for line in commands.read_text().splitlines()]
            for name, path in fixtures.items():
                self.assertEqual(sum(str(path) in args for args in after), 1, name)

    def test_pch_profiles_generation_and_object_reuse(self):
        for pch_enabled in (True, False):
            with self.subTest(pch_enabled=pch_enabled):
                self.check_build(pch_enabled)


if __name__ == "__main__":
    main()

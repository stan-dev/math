"""Check CMake command planning. Run with python3 test/cmake.py."""

import io
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "test" / "prob"))
import getDependencies as get_dependencies
import runTests as runner


def inputs(*flags):
    return runner.processCLIArgs(["test/unit/math/prim/fun", *flags])


class CMakeTest(unittest.TestCase):
    def test_directory_build_with_fresh_filters_and_jobs(self):
        configure, build = runner.build_commands(inputs("--fresh", "-j20", "-f", "sum"),
                                                 ["test/unit/math/prim/fun"], Path("source"), Path("build"))
        self.assertNotIn("--fresh", configure)
        self.assertIn("-DCMAKE_BUILD_TYPE=ReleaseTesting", configure)
        self.assertEqual(build[build.index("--config") + 1], "ReleaseTesting")
        self.assertIn("-DSTAN_TEST_NAME_FILTERS=sum", configure)
        self.assertEqual(build[build.index("--parallel") + 1], "20")
        self.assertEqual(build[build.index("--target") + 1:], ["unit_math_prim_fun_all_tests"])
        execute = runner.test_command("test/unit/math/prim/fun", Path("build"), "Release", {}, False)
        self.assertEqual(Path(execute[0]).name, "unit_math_prim_fun_all_tests")

    def test_fresh_deletes_the_cache_and_keeps_objects(self):
        commands = []

        def run(command, **kwargs):
            commands.append(command)
            return subprocess.CompletedProcess(command, 0, stdout="")

        with tempfile.TemporaryDirectory() as working_dir:
            build = Path(working_dir) / "build"
            objects = build / "test/unit/CMakeFiles/sum.dir"
            for directory in (build / "CMakeFiles/3.16.3", objects):
                directory.mkdir(parents=True)
            (build / "CMakeCache.txt").write_text("STAN_THREADS:BOOL=ON\n")
            (objects / "sum_test.cpp.o").write_text("")
            with patch.object(runner.subprocess, "run", run), patch.object(sys, "stdout", io.StringIO()):
                code = runner.main(["--configure-only", "--fresh"], ROOT, Path(working_dir), False, io.StringIO())
            self.assertEqual(code, 0)
            self.assertFalse((build / "CMakeCache.txt").exists())
            self.assertFalse((build / "CMakeFiles").exists())
            self.assertTrue((objects / "sum_test.cpp.o").exists())
        self.assertEqual(len(commands), 1)
        self.assertNotIn("--fresh", commands[0])

    def test_build_only_and_test_only(self):
        build = runner.build_commands(inputs("--make-only"), ["test/unit/math/prim/fun"],
                                      Path("source"), Path("build"))
        self.assertEqual(len(build), 2)
        self.assertNotIn("--fresh", build[0])
        self.assertEqual(runner.build_commands(inputs("--test-only"), ["test/unit/math/prim/fun"],
                                              Path("source"), Path("build")), [])

    def test_cmake_cache_controls_configuration_and_mpi_launch(self):
        cache = runner.parse_cache("// comment\nCMAKE_CONFIGURATION_TYPES:STRING=Debug;Release\n"
                                   "STAN_EFFECTIVE_MPI:INTERNAL=ON\nMPIEXEC_PREFLAGS:STRING=--oversubscribe\n")
        command = runner.test_command("test/unit/math/prim/fun", Path("build"), "Debug", cache, True)
        self.assertEqual(command[:4], ["mpirun", "-np", "2", "--oversubscribe"])
        self.assertEqual(Path(command[4]), Path("build/test/unit/Debug/unit_math_prim_fun_all_tests.exe"))

    def test_cli_reads_explicit_arguments(self):
        self.assertEqual(inputs().j, 1)
        self.assertEqual(inputs("-j20").j, 20)
        parsed = runner.processCLIArgs(["build-tests", "build", "--source-regex", "^test/unit/", "--jobs", "4"])
        self.assertEqual((parsed.command, parsed.jobs), ("build-tests", 4))
        self.assertEqual(inputs("--test-jobs", "6").test_jobs, 6)
        for flags in (("-j0",), ("--test-jobs", "0"), ("--fresh", "--test-only"),
                      ("--make-only", "--test-only"), ("--configure-only", "--test-only")):
            with self.subTest(flags=flags), patch.object(sys, "stderr"):
                with self.assertRaises(SystemExit) as error:
                    inputs(*flags)
                self.assertEqual(error.exception.code, 2)

    def test_parallel_test_runs_print_in_order_and_report_failures(self):
        codes = {"unit_math_prim_fun_all_tests": 1, "unit_math_rev_fun_all_tests": 0}

        def run(command, **kwargs):
            name = Path(command[0]).name
            return subprocess.CompletedProcess(command, codes[name], stdout=name + " output\n")

        output = io.StringIO()
        with tempfile.TemporaryDirectory() as working_dir, \
                patch.object(runner.subprocess, "run", run), patch.object(sys, "stdout", output):
            code = runner.main(["--test-only", "--test-jobs", "2", "test/unit/math/prim/fun",
                                "test/unit/math/rev/fun"], ROOT, Path(working_dir), False, io.StringIO())
        self.assertEqual(code, 1)
        printed = output.getvalue()
        self.assertLess(printed.index("unit_math_prim_fun_all_tests output"),
                        printed.index("unit_math_rev_fun_all_tests output"))

    def test_trailing_slash_selects_distribution_directories(self):
        commands = []

        def run(command, **kwargs):
            commands.append(command)
            return subprocess.CompletedProcess(command, 0, stdout="")

        with tempfile.TemporaryDirectory() as working_dir, \
                patch.object(runner.subprocess, "run", run), patch.object(sys, "stdout", io.StringIO()):
            code = runner.main(["--test-only", "--test-jobs", "2", "test/prob/*/"],
                               ROOT, Path(working_dir), False, io.StringIO())
        self.assertEqual(code, 0)
        names = {Path(command[0]).name for command in commands}
        self.assertIn("prob_math_bernoulli_all_tests", names)
        self.assertEqual(len(names), len([path for path in (ROOT / "test/prob").iterdir() if path.is_dir()]))

    def test_selection_removes_overlap_and_rejects_invalid_paths(self):
        directory = "test/unit/math/prim/fun"
        matches = {directory: [directory], "sum": [directory + "/sum_test.cpp"]}
        self.assertEqual(runner.select_paths(matches, {directory}), [directory])
        self.assertEqual(runner.target_for(directory + "/sum_test.cpp")[0],
                         "test_unit_math_prim_fun_sum_test")
        with self.assertRaises(ValueError):
            runner.select_paths({"stan/math": ["stan/math"]}, {"stan/math"})
        with self.assertRaises(ValueError):
            runner.select_paths({"missing": []}, set())

    def test_changed_headers_select_existing_dependent_tests(self):
        existing = {"test/unit/math/" + layer + "/fun/sum_test.cpp" for layer in ("prim", "rev", "fwd", "mix")}
        selected = runner.findChangedTests(["stan/math/prim/fun/sum.hpp"], existing)
        self.assertEqual(set(selected), existing)
        self.assertEqual(runner.findChangedTests(["stan/math/prim/fun/sum.hpp"], existing, limit=1), [])
        self.assertRegex(selected[0], runner.source_filter(selected))

    def test_ci_build_reuses_registered_cmake_labels(self):
        inventory = {"tests": [{"name": "sum", "properties": [{"name": "LABELS", "value":
                     ["prim", "test/unit/math/prim/fun/sum_test.cpp"]}]}, {"name": "grad"}]}
        self.assertEqual(runner.build_tests(inventory, "^test/unit/math/prim/", "build", 4, 500),
                         [["cmake", "--build", "build", "--parallel", "4", "--target", "sum"]])
        with self.assertRaisesRegex(ValueError, "No registered tests match"):
            runner.build_tests({"tests": []}, "^test/unit/", "build", 1, 500)

    def test_files_and_fixtures_build_through_one_selected_target(self):
        fixture = "test/prob/bernoulli/bernoulli_test.hpp"
        source = "test/unit/math/prim/fun/sum_test.cpp"
        configure, build = runner.build_commands(inputs(), ["test/unit/math/rev", source, fixture],
                                                 Path("source"), Path("build"))
        self.assertIn("-DSTAN_BUILD_PROB_TESTS=ON", configure)
        self.assertIn("-DCMAKE_EXPORT_COMPILE_COMMANDS=ON", configure)
        self.assertEqual(build[build.index("--target") + 1:],
                         ["unit_math_rev_all_tests", "stan_selected_tests"])
        regex = next(argument for argument in configure
                     if argument.startswith("-DSTAN_SELECTED_TEST_REGEX="))[len("-DSTAN_SELECTED_TEST_REGEX="):]
        self.assertRegex(source, regex)
        self.assertIsNone(re.search(regex, "test/unit/math/prim/fun/cumulative_sum_test.cpp"))
        for mode in ("v", "ffv", "vv", "fd", "fv", "ffd"):
            self.assertRegex(f"test/prob/bernoulli/bernoulli_00012_generated_{mode}_test.cpp", regex)
        self.assertIsNone(re.search(regex, "test/prob/bernoulli/bernoulli_cdf_00000_generated_v_test.cpp"))
        directory_only = runner.build_commands(inputs(), ["test/unit/math/rev"], Path("source"), Path("build"))
        self.assertIn("-DSTAN_SELECTED_TEST_REGEX=", directory_only[0])
        self.assertNotIn("stan_selected_tests", directory_only[1])
        only_configure = runner.build_commands(inputs("--configure-only"), [], Path("source"), Path("build"))
        self.assertEqual(len(only_configure), 1)
        self.assertEqual(only_configure[0][:2], ["cmake", "-S"])

    def test_fixtures_run_their_generated_tests(self):
        labels = ["test/prob/bernoulli/bernoulli_00000_generated_v_test.cpp",
                  "test/prob/bernoulli/bernoulli_00001_generated_fd_test.cpp",
                  "test/prob/bernoulli/bernoulli_cdf_00000_generated_v_test.cpp"]
        inventory = {"tests": [{"name": label, "properties": [{"name": "LABELS", "value":
                     ["prob_math", "bernoulli", label]}]} for label in labels]}
        self.assertEqual(runner.fixture_tests(inventory, ["test/prob/bernoulli/bernoulli_test.hpp"]),
                         labels[:2])
        commands = []

        def run(command, **kwargs):
            commands.append(command)
            if command[0] == "ctest":
                return subprocess.CompletedProcess(command, 0, stdout=json.dumps(inventory))
            return subprocess.CompletedProcess(command, 0, stdout="")

        with tempfile.TemporaryDirectory() as working_dir, \
                patch.object(runner.subprocess, "run", run), patch.object(sys, "stdout", io.StringIO()):
            code = runner.main(["--test-only", "--test-jobs", "2", "test/prob/bernoulli/bernoulli_test.hpp"],
                               ROOT, Path(working_dir), False, io.StringIO())
        self.assertEqual(code, 0)
        self.assertEqual({Path(command[0]).name for command in commands[1:]},
                         {"test_prob_bernoulli_bernoulli_00000_generated_v_test",
                          "test_prob_bernoulli_bernoulli_00001_generated_fd_test"})

    def test_distribution_dependencies_use_the_build_compile_command(self):
        compile_commands = [{"file": "/src/lib/gtest_main.cc", "command":
                             "clang++ -DSTAN -I/src -O0 -MD -MT x.o -MF x.d -o x.o -c /src/lib/gtest_main.cc"}]
        command = get_dependencies.dependency_command(compile_commands)
        self.assertEqual(command, ["clang++", "-DSTAN", "-I/src", "-O0", "-MM", "-MG", "-x", "c++"])
        output = ("normal_test.o: test/prob/normal/normal_test.hpp \\\n"
                  "  /src/stan/math/prim/prob/normal_lpdf.hpp /usr/include/eigen.h \\\n"
                  "  /src/stan/math/prim/fun/constants.hpp\n")
        self.assertEqual(get_dependencies.parse_dependencies(output, Path("/src")),
                         {"stan/math/prim/prob/normal_lpdf.hpp", "stan/math/prim/fun/constants.hpp"})

    def test_ci_build_aggregates_only_complete_selected_groups(self):
        inventory = {"tests": [{"name": name, "properties": [{"name": "LABELS", "value":
                     ["test/unit/" + name + ".cpp", "stan_build:all"]
                     + (["stan_build:pair"] if name != "three" else [])}]} for name in ("one", "two", "three")]}
        command = ["cmake", "--build", "build", "--parallel", "4", "--target"]
        self.assertEqual(runner.build_tests(inventory, "^test/unit/", "build", 4, 500), [command + ["all"]])
        self.assertEqual(runner.build_tests(inventory, "^test/unit/(one|two)", "build", 4, 500), [command + ["pair"]])
        self.assertEqual(runner.build_tests(inventory, "one", "build", 4, 500), [command + ["one"]])


if __name__ == "__main__":
    unittest.main()

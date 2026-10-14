"""Regression checks for directory grouping and runner failure handling."""

from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import cmake_runner

sys.path.insert(0, str(cmake_runner.ROOT))
import runTests


class CMakeRunnerTest(unittest.TestCase):
    def inputs(self, build_dir, **overrides):
        values = dict(tests=["./test/unit/math/prim/fun"], changed=False,
                      debug=False, f=[], only_functions=[], build_dir=build_dir,
                      test_only=False, make_only=False, config="Release", e=-1,
                      j=6, cmake_arg=[], run_all=False)
        values.update(overrides)
        return SimpleNamespace(**values)

    def test_directory_builds_and_runs_one_group(self):
        with tempfile.TemporaryDirectory() as build_dir, patch.object(cmake_runner, "command") as command:
            self.assertEqual(cmake_runner.run(self.inputs(build_dir)), 0)
            configure, build, execute = [call.args[0] for call in command.call_args_list]
            self.assertIn("-DSTAN_TEST_FILTER=", configure)
            self.assertEqual(build[build.index("--target") + 1:], ["unit_math_prim_fun_all_tests"])
            self.assertEqual(Path(execute[0]).name, "unit_math_prim_fun_all_tests")
            self.assertIn("test-results/unit_math_prim_fun_all_tests.xml", execute[1])

    def test_file_builds_individual_executable(self):
        self.assertEqual(cmake_runner.target_for("test/unit/math/prim/fun/sum_test.cpp"),
                         ("test_unit_math_prim_fun_sum_test", "test/unit"))

    def test_other_directory_targets(self):
        for path, target in [("test/unit", "unit_all_tests"),
                             ("test/unit/math", "unit_math_all_tests"),
                             ("test/prob/bernoulli", "prob_math_bernoulli_all_tests"),
                             ("test/expressions", "expressions_all_tests")]:
            self.assertEqual(cmake_runner.target_for(path)[0], target)

    def test_overlap_runs_only_parent(self):
        self.assertEqual(cmake_runner.select_paths([
            "test/unit/math/prim/fun/sum_test.cpp", "./test/unit/math/prim/fun",
            "test/unit/math/prim", "test/unit/math/prim"]), ["test/unit/math/prim"])

    def test_filters_do_not_turn_directory_into_individual_targets(self):
        with tempfile.TemporaryDirectory() as build_dir, patch.object(cmake_runner, "command") as command:
            cmake_runner.run(self.inputs(build_dir, f=["sum", "prim"], make_only=True))
            self.assertIn("-DSTAN_TEST_NAME_FILTERS=sum;prim", command.call_args_list[0].args[0])
            self.assertEqual(command.call_args_list[1].args[0][-1], "unit_math_prim_fun_all_tests")
            self.assertEqual(command.call_count, 2)

    def test_test_only_skips_build(self):
        with tempfile.TemporaryDirectory() as build_dir, patch.object(cmake_runner, "command") as command:
            cmake_runner.run(self.inputs(build_dir, test_only=True))
            self.assertEqual(command.call_count, 1)
            self.assertEqual(Path(command.call_args.args[0][0]).name, "unit_math_prim_fun_all_tests")

    def test_run_all_preserves_failure_and_continues(self):
        with tempfile.TemporaryDirectory() as build_dir, patch.object(cmake_runner, "command") as command:
            command.side_effect = [subprocess.CalledProcessError(1, "failed-test"), None]
            result = cmake_runner.run(self.inputs(build_dir, tests=[
                "test/unit/math/prim/fun/sum_test.cpp",
                "test/unit/math/prim/fun/abs_test.cpp"], test_only=True, run_all=True))
            self.assertEqual(result, 1)
            self.assertEqual(command.call_count, 2)

    def test_build_failure_stops_before_execution(self):
        with tempfile.TemporaryDirectory() as build_dir, patch.object(cmake_runner, "command") as command:
            command.side_effect = [None, subprocess.CalledProcessError(2, "build")]
            with self.assertRaises(subprocess.CalledProcessError):
                cmake_runner.run(self.inputs(build_dir))
            self.assertEqual(command.call_count, 2)

    def test_cli_caps_jobs(self):
        with patch.object(sys, "argv", ["runTests.py", "-j20", "test/unit"]):
            self.assertEqual(runTests.processCLIArgs().j, 6)

    def test_cli_rejects_invalid_job_count(self):
        with patch.object(sys, "argv", ["runTests.py", "-j0", "test/unit"]), patch.object(sys, "stderr"):
            with self.assertRaises(SystemExit) as error:
                runTests.processCLIArgs()
            self.assertEqual(error.exception.code, 2)

    def test_legacy_backend_dispatch(self):
        with patch.object(sys, "argv", ["runTests.py", "--backend", "make", "test/unit"]), patch.object(runTests, "runMake") as run:
            runTests.main()
            self.assertEqual(run.call_args.args[0].tests, ["test/unit"])

    def test_invalid_path_rejected(self):
        with self.assertRaises(ValueError):
            cmake_runner.select_paths(["stan/math"])

    def test_filtered_file_selection_cannot_succeed_without_tests(self):
        with self.assertRaises(ValueError):
            cmake_runner.run(self.inputs("build", tests=[
                "test/unit/math/prim/fun/sum_test.cpp"], f=["cholesky"], test_only=True))

    def test_multi_configuration_executable_path(self):
        with tempfile.TemporaryDirectory() as build_dir, patch.object(cmake_runner, "command") as command:
            Path(build_dir, "CMakeCache.txt").write_text("CMAKE_CONFIGURATION_TYPES:STRING=Debug;Release\n")
            cmake_runner.run(self.inputs(build_dir, test_only=True, config="Debug"))
            self.assertEqual(Path(command.call_args.args[0][0]),
                             Path(build_dir, "test/unit/Debug/unit_math_prim_fun_all_tests"))


if __name__ == "__main__":
    unittest.main()

"""Check CMake command planning. Run with python3 test/cmake.py."""

from pathlib import Path
import re
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import runTests as runner


def inputs(*flags):
    return runner.processCLIArgs(["test/unit/math/prim/fun", *flags])


class CMakeTest(unittest.TestCase):
    def test_directory_build_with_fresh_filters_and_jobs(self):
        configure, build = runner.build_commands(inputs("--fresh", "-j20", "-f", "sum"),
                                                 ["test/unit/math/prim/fun"], Path("source"), Path("build"))
        self.assertIn("--fresh", configure)
        self.assertIn("-DSTAN_TEST_NAME_FILTERS=sum", configure)
        self.assertEqual(build[build.index("--parallel") + 1], "20")
        self.assertEqual(build[build.index("--target") + 1:], ["unit_math_prim_fun_all_tests"])
        execute = runner.test_command("test/unit/math/prim/fun", Path("build"), "Release", {}, False)
        self.assertEqual(Path(execute[0]).name, "unit_math_prim_fun_all_tests")

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
        for flags in (("-j0",), ("--fresh", "--test-only"), ("--make-only", "--test-only")):
            with self.subTest(flags=flags), patch.object(sys, "stderr"):
                with self.assertRaises(SystemExit) as error:
                    inputs(*flags)
                self.assertEqual(error.exception.code, 2)

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

    def test_distribution_selection(self):
        self.assertEqual(runner.changed_distributions(["test/unit/math/prim/fun/sum_test.cpp"]), "")
        self.assertEqual(runner.changed_distributions(["stan/math/prim/fun/sum.hpp"]), "^test/prob/")
        selected = runner.changed_distributions(["test/prob/bernoulli/bernoulli_test.hpp"])
        for mode in ("v", "ffv", "vv", "fd", "fv", "ffd"):
            self.assertRegex(f"test/prob/bernoulli/bernoulli_00000_generated_{mode}_test.cpp", selected)
        self.assertIsNone(re.search(selected, "test/prob/bernoulli/bernoulli_cdf_00000_generated_v_test.cpp"))


if __name__ == "__main__":
    unittest.main()

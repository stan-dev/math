import re
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import cmake_ci
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


class CmakeCISelectionTest(unittest.TestCase):
    def test_unit_only_changes_skip_distributions(self):
        self.assertEqual(cmake_ci.changed_distributions(
            ["test/unit/math/prim/fun/sum_test.cpp"]), "")

    def test_fixture_change_selects_all_modes_and_batches(self):
        selected = cmake_ci.changed_distributions(
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
                selected = cmake_ci.changed_distributions([dependency])
                self.assertRegex("test/prob/bernoulli/bernoulli_00000_generated_v_test.cpp", selected)
                self.assertRegex("test/prob/normal/normal_00000_generated_ffv_test.cpp", selected)

    def test_force_all_with_no_changed_files(self):
        self.assertTrue(cmake_ci.changed_distributions([], run_all=True))

    def test_changed_units_preserve_runner_selection(self):
        with patch("runTests.findChangedTests", return_value=[
            "test/unit/math/prim/core/get_num_threads_test"]):
            selected = cmake_ci.changed_units()
        self.assertRegex("test/unit/math/prim/core/get_num_threads_test.cpp", selected)
        self.assertIsNone(re.search(selected, "test/unit/math/prim/core/set_num_threads_test.cpp"))


if __name__ == "__main__":
    unittest.main()

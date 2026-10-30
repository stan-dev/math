import varmat_compatibility
import json
import tempfile
import unittest
import os
import threading
from unittest.mock import patch


class VarmatCompatibilityMainTest(unittest.TestCase):
    def test_initial_build_failure_propagates_compiler_output(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(
                varmat_compatibility, "WORKING_FOLDER", directory), patch.object(
                varmat_compatibility, "_compile_worker", threading.local()), patch.object(
                varmat_compatibility, "run_command", side_effect=[
                    (True, b"configured", b""),
                    (False, b"PCH build failed", b"compiler diagnostic")]) as command:
            with self.assertRaisesRegex(RuntimeError, "PCH build failedcompiler diagnostic"):
                varmat_compatibility.compile_signature("unused.cpp")
            self.assertEqual(command.call_count, 2)
            self.assertEqual(command.call_args.args[0][-1], "stan_compile_test")

    def test_configuration_failure_propagates(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(
                varmat_compatibility, "get_signatures", return_value=["sum(matrix) => real\n"]), patch.object(
                varmat_compatibility, "build_signature", side_effect=RuntimeError("configuration failed")):
            with self.assertRaisesRegex(RuntimeError, "CMake compilation setup failed"):
                varmat_compatibility.main(["sum"], os.path.join(directory, "results.json"), 1, False)

    def test_sum_works(self):
        f = tempfile.NamedTemporaryFile("w", suffix="_test.json", delete=False)
        f.close()

        varmat_compatibility.main(["sum"], f.name, 1, False)

        with open(f.name) as f:
            results = json.load(f)

        self.assertSetEqual(
            set(
                [
                    "compatible_signatures",
                    "incompatible_signatures",
                    "irrelevant_signatures",
                ]
            ),
            set(results.keys()),
        )
        self.assertSetEqual(
            set(
                [
                    "sum(row_vector) => real\n",
                    "sum(vector) => real\n",
                    "sum(matrix) => real\n",
                ]
            ),
            set(results["compatible_signatures"]),
        )
        self.assertEqual(len(results["incompatible_signatures"]), 0)
        self.assertSetEqual(
            set(
                [
                    "sum(array[] int) => int\n",
                    "sum(array[] real) => real\n",
                    "sum(complex_row_vector) => complex\n",
                    "sum(array[] complex) => complex\n",
                    "sum(complex_matrix) => complex\n",
                    "sum(complex_vector) => complex\n",
                ]
            ),
            set(results["irrelevant_signatures"]),
        )

        os.remove(f.name)


if __name__ == "__main__":
    unittest.main()

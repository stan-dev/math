import varmat_compatibility
import json
import tempfile
import unittest
import os

class VarmatCompatibilityMainTest(unittest.TestCase):
    def test_sum_works(self):
        f = tempfile.NamedTemporaryFile("w", suffix = "_test.json", delete = False)
        f.close()

        varmat_compatibility.main(["sum"], f.name, 1, False)

        with open(f.name) as f:
            results = json.load(f)

        self.assertSetEqual(set(["compatible_signatures", "incompatible_signatures", "irrelevant_signatures"]), set(results.keys()))
        self.assertSetEqual(set(["sum(row_vector) => real\n", "sum(vector) => real\n", "sum(matrix) => real\n"]), set(results["compatible_signatures"]))
        self.assertEqual(len(results["incompatible_signatures"]), 0)
        self.assertSetEqual(set(["sum(array[] int) => int\n", "sum(array[] real) => real\n"]), set(results["irrelevant_signatures"]))

        os.remove(f.name)

if __name__ == '__main__':
    unittest.main()

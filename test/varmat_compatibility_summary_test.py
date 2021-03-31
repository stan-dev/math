from varmat_compatibility_summary import convert_signatures_list_to_functions, select_signatures_matching_functions, remove_signatures_matching_functions, process_results
import unittest

class HelpersTest(unittest.TestCase):
    def setUp(self):
        self.signatures_list = ["chicken(array[] real) => real", "squirrel(array[] real) => real", "frog(array[] real) => real"]
    
    def test_convert(self):
        self.assertSetEqual(set(["chicken", "squirrel", "frog"]), convert_signatures_list_to_functions(self.signatures_list))
    
    def test_select(self):
        self.assertSetEqual(set(["chicken(array[] real) => real"]), select_signatures_matching_functions(self.signatures_list, ["chicken"]))
    
    def test_remove(self):
        self.assertSetEqual(set(["squirrel(array[] real) => real", "frog(array[] real) => real"]), remove_signatures_matching_functions(self.signatures_list, ["chicken"]))

class ProcessResultsTest(unittest.TestCase):
    def setUp(self):
        self.results = {
            "compatible_signatures" : ["chicken(matrix) => real", "dog(matrix) => real"],
            "incompatible_signatures" : ["squirrel(vector) => real", "dog(vector) => real"],
            "irrelevant_signatures" : ["chicken(array[] real) => real", "squirrel(array[] real) => real", "frog(array[] real) => real"]
        }

    def test_which(self):
        self.assertSetEqual(set(["chicken(matrix) => real", "dog(matrix) => real"]), process_results(self.results, functions = [], which = "compatible", fully = False, names = False))
        self.assertSetEqual(set(["squirrel(vector) => real", "dog(vector) => real"]), process_results(self.results, functions = [], which = "incompatible", fully = False, names = False))
        self.assertSetEqual(set(["chicken(array[] real) => real", "squirrel(array[] real) => real", "frog(array[] real) => real"]), process_results(self.results, functions = [], which = "irrelevant", fully = False, names = False))

    def test_fully(self):
        self.assertSetEqual(set(["chicken(matrix) => real"]), process_results(self.results, functions = [], which = "compatible", fully = True, names = False))
        self.assertSetEqual(set(["squirrel(vector) => real"]), process_results(self.results, functions = [], which = "incompatible", fully = True, names = False))
        self.assertSetEqual(set(["frog(array[] real) => real"]), process_results(self.results, functions = [], which = "irrelevant", fully = True, names = False))

    def test_names(self):
        self.assertSetEqual(set(["chicken"]), process_results(self.results, functions = [], which = "compatible", fully = True, names = True))
        self.assertSetEqual(set(["squirrel"]), process_results(self.results, functions = [], which = "incompatible", fully = True, names = True))
        self.assertSetEqual(set(["frog"]), process_results(self.results, functions = [], which = "irrelevant", fully = True, names = True))

    def test_functions(self):
        self.assertSetEqual(set(["chicken(matrix) => real"]), process_results(self.results, functions = ["chicken"], which = "compatible", fully = False, names = False))
        self.assertSetEqual(set(["squirrel(vector) => real"]), process_results(self.results, functions = ["squirrel"], which = "incompatible", fully = False, names = False))
        self.assertSetEqual(set(["frog(array[] real) => real"]), process_results(self.results, functions = ["frog"], which = "irrelevant", fully = False, names = False))

    def test_functions_names(self):
        self.assertSetEqual(set(["chicken"]), process_results(self.results, functions = ["chicken"], which = "compatible", fully = False, names = True))
        self.assertSetEqual(set(["squirrel"]), process_results(self.results, functions = ["squirrel"], which = "incompatible", fully = False, names = True))
        self.assertSetEqual(set(["frog"]), process_results(self.results, functions = ["frog"], which = "irrelevant", fully = False, names = True))

if __name__ == '__main__':
    unittest.main()
import unittest
import subprocess
import sys

error_template = """
Expected "{}" in stdout!

actual stdout:
{}


actual stderr:
{}
"""


class TestExpressionTestingFramework(unittest.TestCase):
    def runCommand(self, command):
        p1 = subprocess.Popen(command, stdout = subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        out, err = p1.communicate()
        return p1.returncode, out, err

    def assertStdoutContains(self, content, stdout, stderr):
        self.assertTrue(content in stdout, msg = error_template.format(content, stdout, stderr))

    def testExpressionNotAcceptedFailure(self):
        return_code, stdout, stderr = self.runCommand((sys.executable, "./runTests.py", "./test/expressions", "-e1", "--make-only", "--only-functions", "bad_no_expressions(matrix) => real"))
        self.assertNotEqual(return_code, 0)

    def testMultipleEvaluationsFailure(self):
        return_code, stdout, stderr = self.runCommand((sys.executable, "./runTests.py", "./test/expressions", "-e1", "--only-functions", "bad_multiple_evaluations(matrix) => matrix"))
        self.assertNotEqual(return_code, 0)
        self.assertStdoutContains("[  FAILED  ] ExpressionTestPrim.bad_multiple_evaluations0", stdout, stderr)
        self.assertStdoutContains("[  FAILED  ] ExpressionTestRev.bad_multiple_evaluations0", stdout, stderr)
        self.assertStdoutContains("[  FAILED  ] ExpressionTestFwd.bad_multiple_evaluations0", stdout, stderr)

    def testWrongResultFailure(self):
        return_code, stdout, stderr = self.runCommand((sys.executable, "./runTests.py", "./test/expressions", "-e1", "--only-functions", "bad_wrong_value(matrix) => real"))
        self.assertNotEqual(return_code, 0)
        self.assertStdoutContains("[  FAILED  ] ExpressionTestPrim.bad_wrong_value0", stdout, stderr)
        self.assertStdoutContains("[  FAILED  ] ExpressionTestRev.bad_wrong_value0", stdout, stderr)
        self.assertStdoutContains("[  FAILED  ] ExpressionTestFwd.bad_wrong_value0", stdout, stderr)

    def testWrongDerivativeFailure(self):
        return_code, stdout, stderr = self.runCommand((sys.executable, "./runTests.py", "./test/expressions", "-e1", "--only-functions", "bad_wrong_derivatives(vector) => real"))
        self.assertNotEqual(return_code, 0)
        self.assertStdoutContains("[       OK ] ExpressionTestPrim.bad_wrong_derivatives0", stdout, stderr)
        self.assertStdoutContains("[  FAILED  ] ExpressionTestRev.bad_wrong_derivatives0", stdout, stderr)
        self.assertStdoutContains("[  FAILED  ] ExpressionTestFwd.bad_wrong_derivatives0", stdout, stderr)


if __name__ == '__main__':
    unittest.main()

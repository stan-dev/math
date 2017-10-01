import os
import re
import sys

test_re = re.compile('class\s*(.*)\s*:\s*public\s*(.*)\s*{')
TYPES = {'d': "double",
         'fd': "fvar<double>",
         'ffd': "fvar<fvar<double>>",
         'fv': "fvar<var>",
         'ffv': "fvar<fvar<var>>",
         'v': "var"}

TDSTRING = "typedef boost::mpl::vector<{0}, type_{1}<{2}>> {0}_{3}_{1};"
ITTCSTRING = "INSTANTIATE_TYPED_TEST_CASE_P({0}_{3}_{1}, {2}, {0}_{3}_{1});"

if __name__ == "__main__":
    fn = sys.argv[1]
    with open(fn) as f:
        test_contents = f.read()
    groups = test_re.search(test_contents)
    test_name = groups.group(1).strip()
    fixture_name = groups.group(2).strip()
    if "Test" in fixture_name:
        fixture_name += "Fixture"
    basedir, filename = os.path.split(fn)
    gen_file_name = os.path.join(basedir, filename[:-8] + "generated_test.cpp")
    with open(gen_file_name, "w") as f:
        f.write("#include <test/prob/gen_test_header.hpp>\n")
        f.write("#include <" + fn + ">\n\n")
        for abbr, c_type in TYPES.items():
            for i in range(16):
                f.write(TDSTRING.format(test_name, i, c_type, abbr) + "\n")
                f.write(ITTCSTRING.format(test_name, i, fixture_name, abbr) + "\n")
                f.write("\n")

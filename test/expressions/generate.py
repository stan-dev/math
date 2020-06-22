import urllib.request
import re

# signature_list_location = "https://raw.githubusercontent.com/stan-dev/stanc3/master/test/integration/signatures/stan_math_sigs.expected"
signature_list_location = "./stan_math_sigs.expected"
signature_list_is_url = False
exceptions_list_location = "./stan_math_sigs_exceptions.expected"

args2test = ["matrix", "vector", "row_vector"]
arg_types = {
    'int': "int",
    'int[]': "std::vector<int>",
    'int[,]': "std::vector<std::vector<int>>",
    'real': "SCALAR",
    'real[]': "std::vector<SCALAR>",
    'real[,]': "std::vector<std::vector<SCALAR>>",
    'vector': "Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>",
    'vector[]': "std::vector<Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>>",
    'row_vector': "Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>",
    'row_vector[]': "std::vector<Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>>",
    'matrix': "Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>",
    '(vector, vector, data real[], data int[]) => vector': "auto",
}

def parse_signature(signature):
    return_type, rest = signature.split(" ", 1)
    function_name, rest = rest.split("(", 1)
    args = re.findall(r"(?:[(][^()]+[)][^,()]+)|(?:[^,()]+(?:,*[]])?)", rest)
    args = [i.strip() for i in args if i.strip()]
    return return_type, function_name, args

def make_arg_code(arg, scalar, var_name):
    arg_type = arg_types[arg].replace("SCALAR", scalar)
    if arg == '(vector, vector, data real[], data int[]) => vector':
        return "  %s %s = [](const auto& a, const auto&, const auto&, const auto&){return a;};" % (arg_type,var_name)
    else:
        return "  %s %s = stan::test::make_arg<%s>();\n" % (arg_type, var_name, arg_type)


part_sig = ""
ignored = set()
for signature in open(exceptions_list_location):
    if not signature.endswith(")\n"):
        part_sig += signature
        continue
    part_sig = ""
    ignored.add(signature)

if signature_list_is_url:
    signatures = urllib.request.urlopen(signature_list_location)
else:
    signatures = open(signature_list_location)

test_n = {}
tests=[]
for signature in signatures:
    if not signature.endswith(")\n"):
        part_sig += signature
        continue
    return_type, function_name, args = parse_signature(part_sig + signature)
    part_sig = ""
    if signature in ignored:
        continue
    for arg2test in args2test:
        if arg2test in args:
            break
    else:
        continue
    func_test_n = test_n.get(function_name,0)
    test_n[function_name] = func_test_n + 1

    test_code = ""
    for overload, scalar in (("Prim", "double"), ("Rev", "stan::math::var"), ("Fwd", "stan::math::fvar<double>")):
        test_code += "TEST(ExpressionTest%s, %s%d){\n" % (overload, function_name, func_test_n)

        for n, arg in enumerate(args):
            test_code += make_arg_code(arg, scalar, "arg_mat%d" % n)
        test_code += "  auto res_mat = stan::math::%s(" % function_name
        for n in range(len(args) - 1):
            test_code += "arg_mat%d, " % n
        test_code += "arg_mat%d);\n\n" % (len(args) - 1)

        for n, arg in enumerate(args):
            test_code += make_arg_code(arg, scalar, "arg_expr%d" % n)
        test_code += "  auto res_expr = stan::math::%s(" % function_name
        for n, arg in enumerate(args[:-1]):
            if arg in arg2test:
                test_code += "1*arg_expr%d, " % n
            else:
                test_code += "arg_expr%d, " % n
        if args[-1] in arg2test:
            test_code += "1*arg_expr%d);\n\n" % (len(args) - 1)
        else:
            test_code += "arg_expr%d);\n\n" % (len(args) - 1)

        test_code += "  EXPECT_STAN_EQ(res_expr, res_mat);\n"

        if overload == "Rev":
            test_code += "  (stan::math::sum(res_mat) + stan::math::sum(res_expr)).grad();\n"
            for n, arg in enumerate(args):
                if arg == '(vector, vector, data real[], data int[]) => vector':
                    continue
                test_code += "  EXPECT_STAN_ADJ_EQ(arg_expr%d,arg_mat%d);\n"%(n,n)
        test_code +="}\n\n"
    tests.append(test_code)

with open("tests.cpp","w") as out:
    out.write("#include <test/expressions/expression_test_helpers.hpp>\n\n")
    for test in tests:
        out.write(test)

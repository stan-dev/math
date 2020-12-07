import sys
import subprocess

sys.path.append("test")
from sig_utils import *
import itertools

build_folder = "./test/benchmark/"

benchmark_template = """
void {benchmark_name}(benchmark::State& state) {{
{setup}
  
  for (auto _ : state) {{
{var_conversions}
    
    auto start = std::chrono::high_resolution_clock::now();
    
{code}
    
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds =
      std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    state.SetIterationTime(elapsed_seconds.count());
    stan::math::recover_memory();
  }}
}}
BENCHMARK({benchmark_name})->RangeMultiplier({multi})->Range(1, {max_size})->UseManualTime();

"""

overload_scalar = {
    "Prim": "double",
    "Rev": "stan::math::var",
    "Fwd": "stan::math::fvar<double>",
    "Mix": "stan::math::fvar<stan::math::var>",
}

def run(command):
    print()
    print(command)
    p1 = subprocess.Popen(command)
    p1.wait()

def build(exe_filepath):
    command = make + " " + exe_filepath
    run(command)

def main(functions_or_sigs, cpp_filename="benchmark.cpp", overloads=("Prim", "Rev"), multiplier_param=None,
         max_size_param=None, max_dimm=3):
    all_signatures = get_signatures()
    functions, signatures = handle_function_list(functions_or_sigs)
    functions = set(functions)
    signatures = set(signatures)
    remaining_functions = set(functions)
    parsed_signatures = []
    for signature in all_signatures:
        return_type, function_name, stan_args = parse_signature(signature)
        if signature in signatures or function_name in functions:
            parsed_signatures.append([return_type, function_name, stan_args])
            remaining_functions.discard(function_name)
    if remaining_functions:
        raise NameError("Functions not found: " + ", ".join(remaining_functions))

    result = ""
    for return_type, function_name, stan_args in parsed_signatures:
        dimm = 0
        for arg in stan_args:
            arg_dimm = 0
            if "vector" in arg:
                arg_dimm = 1
            if "matrix" in arg:
                arg_dimm = 2
            if "[" in arg:
                arg_dimm += len(arg.split("[")[1])
            dimm = max(dimm, arg_dimm)
        if dimm > max_dimm:
            continue
        if max_size_param is None:
            if dimm == 0:  # signature with only scalar arguments
                max_size = 1
            else:
                max_size = 1024 * 1024 * 16
                max_size = int(max_size ** (1. / dimm))
        else:
            max_size = max_size_param
        if multiplier_param is None:
            multiplier = 4
            if dimm >= 2:
                multiplier = 2
        else:
            multiplier = multiplier_param
        cpp_arg_templates = []
        overload_opts = []
        for n, stan_arg in enumerate(stan_args):
            cpp_arg_template = get_cpp_type(stan_arg)
            arg_overload_opts = ["Prim"]
            if "SCALAR" in cpp_arg_template and not (function_name in non_differentiable_args and
                                                     n in non_differentiable_args[function_name]):
                arg_overload_opts = overloads
            cpp_arg_templates.append(cpp_arg_template)
            overload_opts.append(arg_overload_opts)
        for arg_overloads in itertools.product(*overload_opts):
            # generate one benchmark
            benchmark_name = function_name
            setup = ""
            var_conversions = ""
            code = "    auto res = stan::math::{}(".format(function_name)
            for n, (arg_overload, cpp_arg_template, stan_arg), in enumerate(
                    zip(arg_overloads, cpp_arg_templates, stan_args)):

                if stan_arg.endswith("]"):
                    stan_arg2, vec = stan_arg.split("[")
                    benchmark_name += "_" + arg_overload + stan_arg2 + str(len(vec))
                else:
                    benchmark_name += "_" + arg_overload + stan_arg
                scalar = overload_scalar[arg_overload]
                arg_type = cpp_arg_template.replace("SCALAR", scalar)
                var_name = "arg" + str(n)
                if function_name in special_arg_values and special_arg_values[function_name][n] is not None:
                    value = special_arg_values[function_name][n]
                else:
                    value = "0.4"
                if scalar == "double":
                    setup += "  {} {} = stan::test::make_arg<{}>({}, state.range(0));\n".format(
                        arg_type,
                        var_name,
                        arg_type,
                        value,
                    )
                else:
                    var_conversions += "  {} {} = stan::test::make_arg<{}>({}, state.range(0));\n".format(
                        arg_type,
                        var_name,
                        arg_type,
                        value,
                    )
                code += var_name + ", "
            code = code[:-2] + ");\n"
            if "Rev" in arg_overloads:
                code += "    stan::test::recursive_sum(res).grad();\n"
            result += benchmark_template.format(benchmark_name=benchmark_name, setup=setup,
                                                var_conversions=var_conversions, code=code, multi=multiplier,
                                                max_size=max_size)
    cpp_filepath = build_folder + cpp_filename
    with open(cpp_filepath, "w") as o:
        o.write("#include <benchmark/benchmark.h>\n")
        o.write("#include <test/expressions/expression_test_helpers.hpp>\n\n")
        o.write(result)
        o.write("BENCHMARK_MAIN();")
    exe_filepath = cpp_filepath.replace(".cpp", exe_extension)
    build(exe_filepath)
    run(exe_filepath)



if __name__ == "__main__":
    main(["sin"])

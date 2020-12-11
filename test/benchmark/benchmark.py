#!/usr/bin/python
from __future__ import print_function
from argparse import ArgumentParser, RawTextHelpFormatter
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


def run_command(command):
    """
    Runs given command and waits until it finishes executing.
    :param command: command to execute
    """
    print()
    print(command)
    p1 = subprocess.Popen(command)
    if p1.wait() != 0:
        raise RuntimeError("command failed: " + command)


def build(exe_filepath):
    """
    Builds a file using make.
    :param exe_filepath: File to build
    """
    command = make + " " + exe_filepath
    run_command(command)


def run_benchmark(exe_filepath, n_repeats=1, csv_out_file=None):
    """
    Runs a benchmark
    :param exe_filepath: path to the benchmark executable
    :param n_repeats: how many times to repeat each benchmark
    :param csv_out_file: path to csv fle to store benchmark results into
    """
    command = exe_filepath
    if n_repeats > 1:
        command += " --benchmark_repetitions={} --benchmark_report_aggregates_only=true".format(n_repeats)
    if csv_out_file is not None:
        command += " --benchmark_out={} --benchmark_out_format=csv".format(csv_out_file)
    run_command(command)


def plot_results(csv_file, out_file=""):
    """
    Plots benchmark results.
    :param csv_file: path to csv file containing results to plot
    :param out_file: path to image file to store figure into. If empty string opens it in an interactive window.
    """
    import pandas, numpy, matplotlib
    if out_file:
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    data = pandas.read_csv(csv_file)
    data["sig"] = [i.split("/")[0] for i in data["name"]]
    data["size"] = [int(i.split("/")[1]) for i in data["name"]]

    signatures = data["sig"]
    sizes = data["size"]
    times = data["real_time"]

    plt.tight_layout()
    plt.semilogx()
    plt.xlabel("size")
    plt.ylabel("time[us]")

    for signature in numpy.unique(signatures):
        selector = signatures == signature
        sig_sizes = sizes[selector]
        sig_times = times[selector]
        plt.plot(sig_sizes, sig_times, label=signature)

    plt.legend()
    if out_file:
        plt.savefig(out_file)
    else:
        plt.show()


def main(functions_or_sigs, cpp_filename="benchmark.cpp", overloads=("Prim", "Rev"), multiplier_param=None,
         max_size_param=None, max_dim=3, n_repeats=1, csv_out_file=None, plot=False):
    """
    Generates benchmark code, compiles it and runs the benchmark. Optionally plots the results.
    :param functions_or_sigs: List of function names and/or signatures to benchmark
    :param cpp_filename: filename of cpp file to use
    :param overloads: Which overloads to benchmark
    :param multiplier_param: Multiplyer, by which to increase argument size.
    :param max_size_param: Maximum argument size.
    :param max_dim: Maximum number of argument dimensions to benchmark. Signatures with any argument with
             larger number of dimensions are skipped."
    :param n_repeats: Number of times to repeat each benchmark.
    :param csv_out_file: Filename of the csv file to store benchmark results in.
    :param plot: Filename of bmp or csv fle to store plot into. If  filename is empty, opens a window with graph.
    """
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
        if dimm > max_dim:
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
                    var_conversions += "    {} {} = stan::test::make_arg<{}>({}, state.range(0));\n".format(
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
    if plot and csv_out_file is None:
        csv_out_file = ".benchmark.csv"
    run_benchmark(exe_filepath, n_repeats, csv_out_file)
    if plot:
        plot_results(csv_out_file)


def processCLIArgs():
    """
    Define and process the command line interface to the benchmark.py script.
    """
    parser = ArgumentParser(
        description="Generate and run_command benchmarks.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "functions",
        nargs="+",
        type=str,
        default=[],
        help="Signatures and/or function names to benchmark.",
    )
    parser.add_argument(
        "--overloads",
        nargs="+",
        type=str,
        default=["Prim", "Rev"],
        help="Which overload combinations to benchmark. Possible values: Prim, Rev, Fwd, Mix. Defaults to Prim and Rev.",
    )
    parser.add_argument(
        "--multiplier",
        type=int,
        default=None,
        help="Multiplyer, by which to increase argument size. Defaults to 4 for functions with "
             "1-dimensional arguments and 2 for other functions.",
    )
    parser.add_argument(
        "--max_size",
        type=int,
        default=None,
        help="Maximum argument size. Defaults to (16000000)**(1/dimm), where dimm is the largest "
             "number of dimensions of arguments."
    )
    parser.add_argument(
        "--max_dim",
        type=int,
        default=3,
        help="Maximum number of argument dimensions to benchmark. Signatures with any argument with "
             "larger number of dimensions are skipped."
    )
    parser.add_argument(
        "--cpp",
        metavar="filename",
        type=str,
        default="benchmark.cpp",
        help="Filename of the cpp file to generate.",
    )
    parser.add_argument(
        "--repeats",
        metavar="N",
        type=int,
        default=1,
        help="Number of times to repeat each benchmark.",
    )
    parser.add_argument(
        "--csv",
        metavar="filename",
        type=str,
        default=None,
        help="Filename of the csv file to store benchmark results in. By default does not store results.",
    )
    parser.add_argument(
        "--plot",
        metavar="filename",
        type=str,
        default=False,
        help="Filename store plotted graph into. If filename is empty, opens a window with the graph."
             " Plotting requires matplotlib and pandas libraries. Default: no plotting.",
    )
    args = parser.parse_args()
    return dict(functions_or_sigs=args.functions, cpp_filename=args.cpp, overloads=args.overloads,
                multiplier_param=args.multiplier, max_size_param=args.max_size, max_dim=args.max_dim,
                csv_out_file=args.csv, n_repeats=args.repeats, plot=args.plot)


if __name__ == "__main__":
    args = processCLIArgs()
    main(**args)

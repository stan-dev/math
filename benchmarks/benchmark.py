#!/usr/bin/python
from __future__ import print_function

import itertools
import numbers
import os
import subprocess
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

HERE = os.path.dirname(os.path.realpath(__file__))
TEST_FOLDER = os.path.abspath(os.path.join(HERE, "..", "test"))
sys.path.append(TEST_FOLDER)

from sig_utils import *

WORKING_FOLDER = "./benchmarks/"

BENCHMARK_TEMPLATE = """
static void {benchmark_name}(benchmark::State& state) {{
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
    benchmark::ClobberMemory();
  }}
}}
BENCHMARK({benchmark_name})->RangeMultiplier({multi})->Range(1, {max_size})->UseManualTime();

"""

CUSTOM_MAIN = """
int main(int argc, char** argv)
{{
  stan::math::ChainableStack::instance_->memalloc_.alloc({});
  stan::math::recover_memory();

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
}}
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
    print(" ".join(command))
    p1 = subprocess.Popen(command)
    if p1.wait() != 0:
        raise RuntimeError("command failed: " + " ".join(command))


def build(exe_filepath):
    """
    Builds a file using make.
    :param exe_filepath: File to build
    """
    run_command([make, exe_filepath])


def run_benchmark(exe_filepath, n_repeats=1, csv_out_file=None):
    """
    Runs a benchmark
    :param exe_filepath: path to the benchmark executable
    :param n_repeats: how many times to repeat each benchmark
    :param csv_out_file: path to csv fle to store benchmark results into
    """
    command = [exe_filepath]
    if n_repeats > 1:
        command.append("--benchmark_repetitions={}".format(n_repeats))
        command.append("--benchmark_display_aggregates_only=true")
    if csv_out_file is not None:
        command.append("--benchmark_out={}".format(csv_out_file))
        command.append("--benchmark_out_format=csv")
    run_command(command)


def pick_color(n):
    str_bit_reversed_n = "{:015b}".format(n + 1)[::-1]
    r = 0.9 * ((int(str_bit_reversed_n[0::3], 2) / 2.0 ** 5 + 0.3) % 1)
    g = 0.9 * ((int(str_bit_reversed_n[1::3], 2) / 2.0 ** 5 + 0.3) % 1)
    b = 0.9 * ((int(str_bit_reversed_n[2::3], 2) / 2.0 ** 5 + 0.3) % 1)
    return r, g, b


def plot_results(csv_filename, out_file="", plot_log_y=False):
    """
    Plots benchmark results.
    :param csv_filename: path to csv file containing results to plot
    :param out_file: path to image file to store figure into. If it equals to "window" opens it in an interactive window.
    """
    import pandas
    import numpy
    import matplotlib

    if out_file != "window":
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with open(csv_filename) as f:
        # google benchmark writes some non-csv data at beginning
        for line in iter(f.readline, ""):
            if line.startswith("name,iterations"):
                f.seek(f.tell() - len(line) - 1, os.SEEK_SET)
                break
        data = pandas.read_csv(f)
    name_split = data["name"].str.split("/", expand=True)
    timing_data = pandas.concat(
        [name_split.iloc[:, :2], data["real_time"]],
        axis=1,
    ).loc[name_split.iloc[:, 2]=="manual_time", :]
    timing_data.columns = ["signatures", "sizes", "times"]
    timing_data.loc[:, "sizes"] = timing_data["sizes"].astype(int)
    timing_data.loc[:, "times"] /= 1000  # convert to microseconds

    fig, ax = plt.subplots(figsize=(10, 10))
    fig.set_tight_layout(True)
    ax.set_xscale("log")
    if plot_log_y:
        ax.set_yscale("log")
    ax.set_xlabel("size")
    ax.set_ylabel("time[us]")
	
    for n, (signature, sub_data) in enumerate(timing_data.groupby("signatures")):
        avg_sig_times = (
            sub_data.groupby(by="sizes")["times"]
                .median()
                .reset_index()
                .sort_values(by="sizes")
        )
        ax.plot(
            avg_sig_times["sizes"],
            avg_sig_times["times"],
            label=signature,
            color=pick_color(n),
        )
    for n, (signature, sub_data) in enumerate(timing_data.groupby("signatures")):
        ax.plot(
            sub_data["sizes"],
            sub_data["times"],
            "x",
            color=pick_color(n),
            label="_nolegend_",
			scaley=False,
        )
    [
        spine.set_visible(False)
        for loc, spine in ax.spines.items()
        if loc in ["top", "right", "left", "bottom"]
    ]
    ax.minorticks_off()
    ax.grid()
    ax.legend()
    if out_file == "window":
        plt.show()
    else:
        fig.savefig(out_file, bbox_inches="tight", dpi=300)


def plot_compare(csv_filename, reference_csv_filename, out_file="", plot_log_y=False):
    """
    Plots benchmark speedup compared to reference results.
    :param csv_filename: path to csv file containing results to plot
    :param reference_csv_filename: path to csv file containing reference results to plot
    :param out_file: path to image file to store figure into. If it equals to "window" opens it in an interactive window.
    """
    import pandas, numpy, matplotlib

    if out_file != "window":
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with open(csv_filename) as f:
        # google benchmark writes some non-csv data at beginning
        for line in iter(f.readline, ""):
            if line.startswith("name,iterations"):
                f.seek(f.tell() - len(line) - 1, os.SEEK_SET)
                break
        data = pandas.read_csv(f)
    with open(reference_csv_filename) as f:
        # google benchmark writes some non-csv data at beginning
        for line in iter(f.readline, ""):
            if line.startswith("name,iterations"):
                f.seek(f.tell() - len(line) - 1, os.SEEK_SET)
                break
        reference_data = pandas.read_csv(f)
    name_split = data["name"].str.split("/", expand=True)
    timing_data = pandas.concat(
        [name_split.iloc[:, :2], data["real_time"]],
        axis=1,
    ).loc[name_split.iloc[:, 2]=="manual_time", :]
    reference_name_split = reference_data["name"].str.split("/", expand=True)
    reference_timing_data = pandas.concat(
        [
            reference_name_split.iloc[:, :2],
            reference_data["real_time"],
        ],
        axis=1,
    ).loc[reference_name_split.iloc[:, 2]=="manual_time", :]
    timing_data.columns = reference_timing_data.columns = [
        "signatures",
        "sizes",
        "times",
    ]

    same_in_last_selector = reference_timing_data["signatures"].isin(
        timing_data["signatures"]
    )
    reference_timing_data = reference_timing_data.loc[same_in_last_selector, :]
    assert (reference_timing_data["signatures"] == timing_data["signatures"]).all()
    assert (reference_timing_data["sizes"] == timing_data["sizes"]).all()

    timing_data["speedup"] = reference_timing_data["times"] / timing_data["times"]
    timing_data["sizes"] = timing_data["sizes"].astype(int)

    fig, ax = plt.subplots(figsize=(10, 10))
    fig.set_tight_layout(True)
    ax.set_xscale("log")
    if plot_log_y:
        ax.set_yscale("log")
    ax.set_xlabel("size")
    ax.set_ylabel("speedup")

    for n, (signature, sub_data) in enumerate(timing_data.groupby("signatures")):
        avg_sig_speedups = (
            sub_data.groupby(by="sizes")["speedup"]
                .median()
                .reset_index()
                .sort_values(by="sizes")
        )
        ax.plot(
            avg_sig_speedups["sizes"],
            avg_sig_speedups["speedup"],
            label=signature,
            color=pick_color(n),
        )
    plt.plot([1, max(timing_data["sizes"])], [1, 1], "--", color="gray")
    for n, (signature, sub_data) in enumerate(timing_data.groupby("signatures")):
        ax.plot(
            sub_data["sizes"],
            sub_data["speedup"],
            "x",
            color=pick_color(n),
            label="_nolegend_",
			scaley=False,
        )

    [
        spine.set_visible(False)
        for loc, spine in ax.spines.items()
        if loc in ["top", "right", "left", "bottom"]
    ]
    ax.minorticks_off()
    ax.grid()
    ax.legend()
    if out_file == "window":
        plt.show()
    else:
        fig.savefig(out_file, bbox_inches="tight", dpi=300)


def benchmark(
        functions_or_sigs,
        cpp_filename="benchmark.cpp",
        overloads=("Prim", "Rev"),
        multiplier_param=None,
        max_size_param=None,
        max_dim=3,
        n_repeats=1,
        skip_similar_signatures=False,
        csv_out_file=None,
        opencl=False,
        varmat=False,
):
    """
    Generates benchmark code, compiles it and runs the benchmark.
    :param functions_or_sigs: List of function names and/or signatures to benchmark
    :param cpp_filename: filename of cpp file to use
    :param overloads: Which overloads to benchmark
    :param multiplier_param: Multiplyer, by which to increase argument size.
    :param max_size_param: Maximum argument size.
    :param max_dim: Maximum number of argument dimensions to benchmark. Signatures with any argument with
             larger number of dimensions are skipped."
    :param n_repeats: Number of times to repeat each benchmark.
    :param skip_similar_signatures: Whether to skip similar signatures. Two signatures are similar if they
             difffer only in similar vector types, which are vector, row_vector and real[].
    :param csv_out_file: Filename of the csv file to store benchmark results in.
    """
    all_signatures = get_signatures()
    functions, signatures = handle_function_list(functions_or_sigs)
    functions = set(functions)
    signatures = set(signatures)
    remaining_functions = set(functions)
    parsed_signatures = []
    ref_signatures = set()
    for signature in all_signatures:
        return_type, function_name, stan_args = parse_signature(signature)
        reference_args = tuple(reference_vector_argument(i) for i in stan_args)
        if (
                skip_similar_signatures
                and (function_name, reference_args) in ref_signatures
        ):
            continue
        if (signature in signatures) or (function_name in functions):
            parsed_signatures.append([return_type, function_name, stan_args])
            remaining_functions.discard(function_name)
            ref_signatures.add((function_name, reference_args))
    for signature in signatures:
        return_type, function_name, stan_args = parse_signature(signature)
        reference_args = tuple(reference_vector_argument(i) for i in stan_args)
        if (
                skip_similar_signatures
                and (function_name, reference_args) in ref_signatures
        ):
            continue
        ref_signatures.add((function_name, reference_args))
        parsed_signatures.append([return_type, function_name, stan_args])
        remaining_functions.discard(function_name)
    if remaining_functions:
        raise NameError(
            "Functions not found: " + ", ".join(sorted(remaining_functions))
        )
    result = ""
    max_args_with_max_dimm = 0
    default_max_size = 1024 * 1024 * 16

    for return_type, function_name, stan_args in parsed_signatures:
        dimm = 0
        args_with_max_dimm = 0
        for arg in stan_args:
            arg_dimm = 0
            if "vector" in arg:
                arg_dimm = 1
            if "matrix" in arg:
                arg_dimm = 2
            if "[" in arg:
                arg_dimm += len(arg.split("]")[0].split("[")[1])
            if arg_dimm == dimm:
                args_with_max_dimm += 1
            elif arg_dimm > dimm:
                dimm = arg_dimm
                args_with_max_dimm = 1
        if dimm > max_dim:
            continue
        max_args_with_max_dimm = max(max_args_with_max_dimm, args_with_max_dimm)
        if max_size_param is None:
            if dimm == 0:  # signature with only scalar arguments
                max_size = 1
            else:
                max_size = default_max_size
                max_size = int(max_size ** (1.0 / dimm))
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
            if "SCALAR" in cpp_arg_template and not (
                    function_name in non_differentiable_args
                    and n in non_differentiable_args[function_name]
            ):
                arg_overload_opts = overloads
            cpp_arg_templates.append(cpp_arg_template)
            overload_opts.append(arg_overload_opts)
        for arg_overloads in itertools.product(*overload_opts):
            # generate one benchmark
            benchmark_name = function_name
            setup = ""
            var_conversions = ""
            if opencl in ("copy", "copy_rev") and return_type not in scalar_stan_types:
                code = "    auto res = stan::math::from_matrix_cl(stan::math::{}(".format(
                    function_name
                )
            else:
                code = "    auto res = stan::math::eval(stan::math::{}(".format(
                    function_name
                )
            for (
                    n,
                    (arg_overload, cpp_arg_template, stan_arg),
            ) in enumerate(zip(arg_overloads, cpp_arg_templates, stan_args)):
                n_vec, inner_type = parse_array(stan_arg)
                if n_vec:
                    benchmark_name += (
                            "_" + arg_overload + "_" + inner_type + str(n_vec)
                    )
                else:
                    benchmark_name += "_" + arg_overload + "_" + stan_arg
                scalar = overload_scalar[arg_overload]
                arg_type = cpp_arg_template.replace("SCALAR", scalar)
                var_name = "arg" + str(n)
                make_arg_function = "make_arg"
                is_argument_autodiff = "var" in arg_type
                is_argument_scalar = stan_arg in scalar_stan_types
                value = 0.4
                if function_name in special_arg_values:
                    if isinstance(special_arg_values[function_name][n], str):
                        make_arg_function = make_special_arg_values[special_arg_values[function_name][n]]
                    elif isinstance(
                            special_arg_values[function_name][n], numbers.Number
                    ):
                        value = special_arg_values[function_name][n]
                if not is_argument_autodiff or (
                      not is_argument_scalar and (
                          opencl == "base" or varmat == "base" or make_arg_function != "make_arg"
                )):
                    arg_type_prim = cpp_arg_template.replace("SCALAR", "double");
                    setup += (
                        "  {} {} = stan::test::{}<{}>({}, state.range(0));\n".format(
                            arg_type_prim,
                            var_name,
                            make_arg_function,
                            arg_type_prim,
                            value,
                        )
                    )
                    if not is_argument_scalar:
                        if opencl == "base" or opencl == "copy_rev":
                            setup += "  auto {} = stan::math::to_matrix_cl({});\n".format(
                                var_name + "_cl", var_name
                            )
                            var_name += "_cl"
                            if is_argument_autodiff:
                                var_conversions += (
                                  "    stan::math::var_value<stan::math::matrix_cl<double>> {}({});\n".format(
                                    var_name + "_var", var_name)
                                )
                                var_name += "_var"
                        elif varmat == "base" and arg_overload == "Rev":
                            var_conversions += "      stan::math::var_value<{}> {}({});\n".format(
                                arg_type_prim, var_name + "_varmat", var_name
                            )
                            var_name += "_varmat"
                        elif is_argument_autodiff: #rev
                            var_conversions += "    {} {} = {};\n".format(
                                arg_type, var_name + "_var", var_name
                            )
                            var_name += "_var"
                else:
                    var_conversions += (
                        "    {} {} = stan::test::{}<{}>({}, state.range(0));\n".format(
                            arg_type,
                            var_name,
                            make_arg_function,
                            arg_type,
                            value,
                        )
                    )
                    if not is_argument_scalar:
                        if opencl == "base" or (opencl == "copy_rev" and not is_argument_autodiff):
                            var_conversions += (
                                "    auto {} = stan::math::to_matrix_cl({});\n".format(
                                    var_name + "_cl", var_name
                                )
                            )
                            var_name += "_cl"
                        elif varmat == "base" and arg_overload == "Rev":
                            var_conversions += (
                                "    auto {} = stan::math::to_var_value({});\n".format(
                                    var_name + "_varmat", var_name
                                )
                            )
                            var_name += "_varmat"
                if (opencl == "copy" or opencl == "copy_rev" and is_argument_autodiff) and not is_argument_scalar:
                    code += "stan::math::to_matrix_cl({}), ".format(var_name)
                elif (
                        varmat == "copy"
                        and not is_argument_scalar
                        and arg_overload == "Rev"
                ):
                    code += "stan::math::to_var_value({}), ".format(var_name)
                else:
                    code += var_name + ", "
            code = code[:-2] + "));\n"
            if "Rev" in arg_overloads:
                code += "    stan::math::grad();\n"
            if opencl == "base":
              code += "    stan::math::opencl_context.queue().finish();\n"
              var_conversions += "    stan::math::opencl_context.queue().finish();\n"
            result += BENCHMARK_TEMPLATE.format(
                benchmark_name=benchmark_name,
                setup=setup,
                var_conversions=var_conversions,
                code=code,
                multi=multiplier,
                max_size=max_size,
            )
    cpp_filepath = os.path.join(WORKING_FOLDER, cpp_filename)
    with open(cpp_filepath, "w") as f:
        f.write("#include <benchmark/benchmark.h>\n")
        f.write("#include <test/expressions/expression_test_helpers.hpp>\n\n")
        f.write(result)
        if "Rev" in overloads:
            # estimate the amount of arena memory the benchmarks will need
            DOUBLE_SIZE = 8
            N_ARRAYS = 4  # vals, adjoints, pointers + 1 for anything else
            f.write(
                CUSTOM_MAIN.format(
                    (max_size_param or default_max_size)
                    * DOUBLE_SIZE
                    * N_ARRAYS
                    * (max_args_with_max_dimm + 1)
                )
            )
        else:
            f.write("BENCHMARK_MAIN();")
    exe_filepath = cpp_filepath.replace(".cpp", exe_extension)
    build(exe_filepath)
    run_benchmark(exe_filepath, n_repeats, csv_out_file)


def main(
        functions_or_sigs,
        cpp_filename="benchmark.cpp",
        overloads=("Prim", "Rev"),
        multiplier_param=None,
        max_size_param=None,
        max_dim=3,
        n_repeats=1,
        skip_similar_signatures=False,
        csv_out_file=None,
        opencl=False,
        varmat=False,
        plot=False,
        plot_log_y=False,
        plot_speedup=False,
        plot_reference=None,
):
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
    :param skip_similar_signatures: Whether to skip similar signatures. Two signatures are similar if they
             difffer only in similar vector types, which are vector, row_vector and real[].
    :param csv_out_file: Filename of the csv file to store benchmark results in.
    :param plot: Filename of bmp or csv fle to store plot into. If  filename is empty, opens a window with graph.
    :param plot_log_y: Use logarithmic y axis for plotting
    :param plot_speedup: plot speedup of OpenCL or varmat overloads compared to CPU ones
    """
    if plot and csv_out_file is None:
        csv_out_file = ".benchmark.csv"
    if plot_speedup and (opencl or varmat):
        if opencl:
            special = "_cl"
        else:
            special = "_varmat"
        opencl_csv_out_file = csv_out_file + special
        if "." in csv_out_file:
            base, ext = csv_out_file.rsplit(".", 1)
            opencl_csv_out_file = base + special + "." + ext
        benchmark(
            functions_or_sigs,
            cpp_filename,
            overloads,
            multiplier_param,
            max_size_param,
            max_dim,
            n_repeats,
            skip_similar_signatures,
            csv_out_file,
            False,
            False,
        )
        benchmark(
            functions_or_sigs,
            cpp_filename,
            overloads,
            multiplier_param,
            max_size_param,
            max_dim,
            n_repeats,
            skip_similar_signatures,
            opencl_csv_out_file,
            opencl,
            varmat,
        )
        plot_compare(opencl_csv_out_file, csv_out_file, plot)
    else:
        benchmark(
            functions_or_sigs,
            cpp_filename,
            overloads,
            multiplier_param,
            max_size_param,
            max_dim,
            n_repeats,
            skip_similar_signatures,
            csv_out_file,
            opencl,
            varmat,
        )
        if plot_reference:
            plot_compare(csv_out_file, plot_reference, plot, plot_log_y)
        elif plot:
            plot_results(csv_out_file, plot, plot_log_y)


class FullErrorMsgParser(ArgumentParser):
    """
    Modified ArgumentParser that prints full error message on any error.
    """

    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)


def processCLIArgs():
    """
    Define and process the command line interface to the benchmark.py script.
    """
    parser = FullErrorMsgParser(
        description="Generate and run_command benchmarks.",
        formatter_class=ArgumentDefaultsHelpFormatter,
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
             "number of dimensions of arguments.",
    )
    parser.add_argument(
        "--max_dim",
        type=int,
        default=3,
        help="Maximum number of argument dimensions to benchmark. Signatures with any argument with "
             "larger number of dimensions are skipped.",
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
        help="Filename store plotted graph into. If filename equals to 'window', opens a window with the graph."
             " Plotting requires matplotlib and pandas libraries. Default: no plotting.",
    )
    parser.add_argument(
        "--plot_log_y",
        default=False,
        action="store_true",
        help="Use logarithmic y axis when plotting.",
    )
    parser.add_argument(
        "--opencl",
        metavar="setting",
        type=str,
        default=False,
        help="Benchmark OpenCL overloads. Possible values: "
             "base - benchmark just the execution time, "
             "copy - include argument copying time"
             "copy_rev - include argument copying time for var arguments only",
    )
    parser.add_argument(
        "--varmat",
        metavar="setting",
        type=str,
        default=False,
        help="Benchmark varmat overloads. Possible values: "
             "base - benchmark just the execution time, "
             "copy - include argument copying time",
    )
    parser.add_argument(
        "--plot_speedup",
        default=False,
        action="store_true",
        help="Plots speedup of OpenCL or varmat overloads compared to Eigen matvar ones. Can only be specified together "
             "with both --plot and either --opencl or --varmat. Cannot be specified together with --plot_reference.",
    )
    parser.add_argument(
        "--plot_reference",
        metavar="filename",
        type=str,
        default=None,
        help="Specify filename of reference run csv output. Plots speedup of this run compared to the reference. "
             "Reference run must have all parameters the same as this one, except possibly --opencl, output files and "
             "plotting parameters. Can only be specified together with --plot. Cannot be specified together with "
             "--plot_cl_speedup.",
    )
    parser.add_argument(
        "--skip_similar_signatures",
        default=False,
        action="store_true",
        help="Skip similar signatures. Two signatures are similar if they"
             "difffer only in similar vector types, which are vector, row_vector and real[].",
    )
    args = parser.parse_args()
    assert not (args.opencl and args.varmat), ValueError(
        "--opencl and --varmat cannot be specified at the same time!"
    )
    if args.plot_reference or args.plot_speedup or args.plot_log_y:
        assert args.plot, ValueError(
            "--plot is required if you specify any of --plot_reference, --plot_speedup, --plot_log_y!"
        )
    main(
        functions_or_sigs=args.functions,
        cpp_filename=args.cpp,
        overloads=args.overloads,
        multiplier_param=args.multiplier,
        max_size_param=args.max_size,
        max_dim=args.max_dim,
        csv_out_file=args.csv,
        n_repeats=args.repeats,
        skip_similar_signatures=args.skip_similar_signatures,
        plot=args.plot,
        plot_log_y=args.plot_log_y,
        opencl=args.opencl,
        plot_speedup=args.plot_speedup,
        plot_reference=args.plot_reference,
        varmat=args.varmat,
    )


if __name__ == "__main__":
    processCLIArgs()

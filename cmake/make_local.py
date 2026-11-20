"""Evaluate supported make/local settings and export CMake defaults."""

import argparse
import os
from pathlib import Path
import re
import shlex
import subprocess


SETTINGS = {
    "CC": "CMAKE_C_COMPILER",
    "CXX": "CMAKE_CXX_COMPILER",
    "O": "STAN_TEST_OPTIMIZATION",
    "OPENCL_DEVICE_ID": "STAN_OPENCL_DEVICE_ID",
    "OPENCL_PLATFORM_ID": "STAN_OPENCL_PLATFORM_ID",
    **{name: name for name in (
        "STAN_THREADS", "STAN_OPENCL", "STAN_MPI", "STAN_NO_RANGE_CHECKS",
        "STAN_TEST_ROW_VECTORS", "STAN_PROB_TEST_ALL")},
    **{name: "STAN_" + name + "_ROOT" for name in (
        "BOOST", "EIGEN", "TBB", "SUNDIALS", "OPENCL", "BENCHMARK")},
    "INTEGRATED_OPENCL": "STAN_INTEGRATED_OPENCL",
}
FLAGS = {
    "CXX": ("CXXFLAGS", "CXXFLAGS_OPTIM", "CXXFLAGS_WARNINGS"),
    "C": ("CFLAGS",),
    "CPP": ("CPPFLAGS", "CPPFLAGS_OPTIM"),
    "LINK": ("LDFLAGS", "LDFLAGS_OPTIM"),
    "LIBRARIES": ("LDLIBS",),
}
FEATURES = {name for name in SETTINGS if name.startswith("STAN_")}


def bracket(value):
    delimiter = "="
    while "]" + delimiter + "]" in value:
        delimiter += "="
    return "[" + delimiter + "[" + value + "]" + delimiter + "]"


def read_local(make, source, local):
    names = list(SETTINGS) + [name for group in FLAGS.values() for name in group]
    escaped_local = str(local).replace("\\", "/").replace(" ", "\\ ").replace("#", "\\#").replace("$", "$$")
    makefile = "include " + escaped_local + "\n"
    for name in names:
        makefile += (f"$(info __STAN_CMAKE_BEGIN_{name}__)\n"
                     f"$(info $(if $(filter file override,$(origin {name})),$({name})))\n"
                     f"$(info __STAN_CMAKE_END_{name}__)\n")
    makefile += ("$(info __STAN_CMAKE_BEGIN_FILES__)\n$(info $(MAKEFILE_LIST))\n"
                 "$(info __STAN_CMAKE_END_FILES__)\n"
                 ".PHONY: stan-cmake-settings\nstan-cmake-settings:\n")
    environment = dict(os.environ)
    for name in ("MAKEFLAGS", "MFLAGS", "GNUMAKEFLAGS", "MAKELEVEL"):
        environment.pop(name, None)
    result = subprocess.run([make, "--no-print-directory", "-s", "-f", "-",
                             "stan-cmake-settings"], cwd=source, env=environment,
                            input=makefile, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    if result.returncode:
        raise ValueError(f"Cannot evaluate {local}:\n{result.stderr}")
    values = {}
    for name in names + ["FILES"]:
        match = re.search(f"__STAN_CMAKE_BEGIN_{name}__\n(.*?)\n__STAN_CMAKE_END_{name}__",
                          result.stdout, re.DOTALL)
        if not match:
            raise ValueError(f"Missing GNU Make output for {name} in {local}")
        values[name] = match.group(1).strip()
    return values


def cmake_settings(values, source):
    output = []
    for name, target in SETTINGS.items():
        value = values[name]
        if name in FEATURES or name == "INTEGRATED_OPENCL":
            value = ("ON" if value else "OFF") if name in FEATURES else ("1" if value else "0")
        elif not value:
            continue
        elif target.endswith("_ROOT"):
            path = Path(value)
            value = str((source / path).resolve())
        elif name in ("CC", "CXX"):
            value = ";".join(shlex.split(value))
        output.append(f"if(NOT DEFINED {target})\n  set({target} {bracket(value)})\nendif()\n")
    for group, names in FLAGS.items():
        value = " ".join(values[name] for name in names if values[name])
        output.append(f"set(STAN_MAKE_LOCAL_{group}_FLAGS {bracket(value)})\n")
    files = []
    words = values["FILES"].split()
    index = 0
    while index < len(words):
        for end in range(index + 1, len(words) + 1):
            path = source / " ".join(words[index:end])
            if path.is_file():
                files.append(str(path.resolve()))
                index = end - 1
                break
        index += 1
    output.append("set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS\n")
    output.extend("  " + bracket(path) + "\n" for path in files)
    output.append(")\n")
    return "".join(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--make", required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--local", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        values = read_local(args.make, args.source, args.local)
        args.output.write_text(cmake_settings(values, args.source))
    except ValueError as error:
        parser.exit(1, str(error) + "\n")


if __name__ == "__main__":
    main()

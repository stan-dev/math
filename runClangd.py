#!/usr/bin/env python3

"""
Generate compile_commands.json for clangd (and the clangd-lsp Claude Code
plugin) with the flags make would use.

The file is gitignored and never committed. Call script with '-h' as an
option to see a helpful message.
"""

from __future__ import annotations

import functools
import json
import os
import re
import shlex
import shutil
import sys
import time
from argparse import ArgumentParser, Namespace, RawTextHelpFormatter
from typing import Any

from utils import ROOT, files_in_folder, run_command, stopErr

CDB_FILE = "compile_commands.json"

# Every path --clean may remove. Nothing outside this list is ever deleted.
GENERATED_PATHS = [
    CDB_FILE,
    CDB_FILE + ".tmp",
    os.path.join(".cache", "clangd"),
]


class Args(Namespace):
    """Parsed command line. cdb options are absent without the cdb command."""

    clean: bool
    command: str | None
    compiler: str
    pin_system_includes: bool
    opencl: bool
    no_tests: bool


def processCLIArgs() -> Args:
    """
    Define and process the command line interface to the runClangd.py script.
    """
    cli_description = (
        "Generate compile_commands.json for clangd. It is gitignored and\n"
        "never committed."
    )
    cli_epilog = (
        "Examples:\n"
        "  ./runClangd.py cdb            write compile_commands.json\n"
        "  ./runClangd.py --clean        remove everything this script made\n"
        "  ./runClangd.py --clean cdb    remove, then regenerate\n"
        "  ./runClangd.py cdb --opencl   OpenCL flags for every entry\n"
    )

    cdb_opts = ArgumentParser(add_help=False)
    cdb_opts.add_argument(
        "--compiler",
        default="clang++",
        help="clang driver written into compile_commands.json (default: clang++)",
    )
    cdb_opts.add_argument(
        "--pin-system-includes",
        action="store_true",
        help="add the system include dirs reported by the compiler as\n"
        "explicit -isystem flags, for when clangd picks up a different\n"
        "standard library than the compiler does",
    )
    cdb_opts.add_argument(
        "--opencl",
        action="store_true",
        help="use the OpenCL flags for every entry, not just stan/math/opencl",
    )
    cdb_opts.add_argument(
        "--no-tests",
        action="store_true",
        help="only write entries for headers, not for test/unit/**.cpp",
    )

    parser = ArgumentParser(
        description=cli_description,
        epilog=cli_epilog,
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="remove every file this script generates, and clangd's index,\n"
        "then run the command if one is given",
    )
    sub = parser.add_subparsers(dest="command")
    sub.add_parser(
        "cdb",
        parents=[cdb_opts],
        help="write compile_commands.json",
        formatter_class=RawTextHelpFormatter,
    )
    args = parser.parse_args(namespace=Args())
    if args.command is None and not args.clean:
        parser.print_help()
        sys.exit(1)
    return args


def git_ok(*args: str) -> bool:
    """True if the git command succeeds."""
    return run_command(["git", *args], capture=True, check=False).returncode == 0


def require_ignored(path: str) -> None:
    """Abort unless git ignores path, so generated files are never committed."""
    if not git_ok("check-ignore", "-q", path):
        stopErr(f"{path} is not gitignored; refusing to generate it", 1)


def clean() -> None:
    """Remove every path in GENERATED_PATHS that exists and git does not track."""
    for rel in GENERATED_PATHS:
        path = os.path.join(ROOT, rel)
        if not os.path.lexists(path):
            continue
        if git_ok("ls-files", "--error-unmatch", rel):
            print(f"refusing to remove {rel}: it is tracked by git")
            continue
        print(f"removing {rel}")
        if os.path.isdir(path) and not os.path.islink(path):
            shutil.rmtree(path)
        else:
            os.remove(path)


def make_vars(names: list[str], opencl: bool = False) -> dict[str, list[str]]:
    """Read make variables through the print-% rule, so make/local is honored."""
    command = ["make", "-s"]
    if opencl:
        command.append("STAN_OPENCL=true")
    command += ["print-" + name for name in names]
    values: dict[str, list[str]] = {}
    for line in run_command(command, capture=True).stdout.splitlines():
        match = re.match(r"^(\w+) = ?(.*)$", line)
        if match and match.group(1) in names:
            values[match.group(1)] = shlex.split(match.group(2))
    missing = [name for name in names if name not in values]
    if missing:
        stopErr(f"make did not print: {', '.join(missing)}", 1)
    return values


@functools.cache
def build_flags(opencl: bool = False) -> tuple[list[str], list[str]]:
    """Compiler flags for Stan headers and for gtest files."""
    v = make_vars(
        ["CXXFLAGS", "CPPFLAGS", "INC_GTEST", "CXXFLAGS_GTEST", "CPPFLAGS_GTEST"],
        opencl=opencl,
    )
    base = v["CXXFLAGS"] + v["CPPFLAGS"]
    gtest = v["INC_GTEST"] + v["CXXFLAGS_GTEST"] + v["CPPFLAGS_GTEST"]
    return base, gtest


def system_includes(compiler: str) -> list[str]:
    """-isystem flags for the dirs the compiler searches for <...> includes."""
    proc = run_command(
        [compiler, "-E", "-v", "-x", "c++", os.devnull], capture=True, check=False
    )
    flags: list[str] = []
    inside = False
    for line in proc.stderr.splitlines():
        if line.startswith("#include <...> search starts here"):
            inside = True
        elif line.startswith("End of search list"):
            break
        elif inside:
            flags += ["-isystem", os.path.normpath(line.strip())]
    return flags


def find_files(top: str, suffix: str) -> list[str]:
    """
    Sorted repo-relative '/' paths under top ending in suffix. top is
    relative to ROOT, which main() makes the working directory.
    """
    return sorted(f.replace(os.sep, "/") for f in files_in_folder(top, suffix))


def is_opencl_path(rel: str) -> bool:
    """True for files that only compile with STAN_OPENCL defined."""
    return rel.startswith(("stan/math/opencl/", "test/unit/math/opencl/"))


def write_cdb(args: Args) -> None:
    """Write compile_commands.json for every header and unit test."""
    require_ignored(CDB_FILE)
    start = time.time()
    pinned = system_includes(args.compiler) if args.pin_system_includes else []

    def entry(rel: str, head: list[str], tail: list[str]) -> dict[str, Any]:
        base, gtest = build_flags(args.opencl or is_opencl_path(rel))
        flags = base + (gtest if tail[0] == "-c" else [])
        return {
            "directory": ROOT,
            "file": os.path.join(ROOT, rel),
            "arguments": [args.compiler, *head, *flags, *pinned, *tail],
        }

    headers = find_files("stan", ".hpp")
    tests = [] if args.no_tests else find_files("test/unit", ".cpp")
    entries = [entry(h, ["-x", "c++-header"], ["-fsyntax-only", h]) for h in headers]
    entries += [entry(t, [], ["-c", t, "-o", t[: -len(".cpp")] + ".o"]) for t in tests]
    tmp = os.path.join(ROOT, CDB_FILE + ".tmp")
    with open(tmp, "w") as f:
        json.dump(entries, f)
    os.replace(tmp, os.path.join(ROOT, CDB_FILE))
    print(
        f"wrote {CDB_FILE}: {len(headers)} headers, {len(tests)} tests"
        f" ({time.time() - start:.1f}s)"
    )
    print("rerun after changing make/local or adding files.")


def main() -> None:
    args = processCLIArgs()
    os.chdir(ROOT)
    if args.clean:
        clean()
    if args.command == "cdb":
        write_cdb(args)


if __name__ == "__main__":
    main()

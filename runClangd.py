#!/usr/bin/env python3

"""
Generate local, gitignored tooling files for clangd and AI coding agents.

  cdb      write compile_commands.json (used by clangd and the clangd-lsp
           Claude Code plugin) with the flags make would use.
  catalog  write an API catalog of Stan Math to .agents/catalog/, one line
           per symbol: name | header | signature | brief
  all      both of the above.

None of the generated files are ever committed. Call script with '-h' as an
option to see a helpful message.
"""

from __future__ import annotations

import bisect
import datetime
import functools
import json
import os
import re
import shlex
import shutil
import sys
import textwrap
import time
from argparse import ArgumentParser, Namespace, RawTextHelpFormatter
from collections import defaultdict
from concurrent.futures import Executor, ProcessPoolExecutor
from dataclasses import dataclass, field
from typing import Any, NamedTuple, TypedDict

from utils import ROOT, files_in_folder, run_command, stopErr

MAX_JOBS = 16
CDB_FILE = "compile_commands.json"
CATALOG_DIR = os.path.join(".agents", "catalog")

# Every path --clean may remove. Nothing outside this list is ever deleted.
GENERATED_PATHS = [
    CDB_FILE,
    CDB_FILE + ".tmp",
    CATALOG_DIR,
    CATALOG_DIR + ".tmp",
    os.path.join(".cache", "clangd"),
]

# libclang cursor kinds recorded in the catalog, by CursorKind name.
KEEP_KINDS = {
    "FUNCTION_TEMPLATE": "function",
    "FUNCTION_DECL": "function",
    "CLASS_TEMPLATE": "class",
    "CLASS_TEMPLATE_PARTIAL_SPECIALIZATION": "class",
    "STRUCT_DECL": "class",
    "CLASS_DECL": "class",
    "TYPE_ALIAS_TEMPLATE_DECL": "alias",
    "TYPE_ALIAS_DECL": "alias",
    "TYPEDEF_DECL": "alias",
    "VAR_DECL": "variable",
}

# The clang python bindings ship without type information, so their module
# and cursor objects are typed as Any.
CIndex = Any
Cursor = Any


class Args(Namespace):
    """Parsed command line. Options of other subcommands are absent."""

    clean: bool
    command: str | None
    compiler: str
    pin_system_includes: bool
    opencl: bool
    no_tests: bool
    j: int
    out: str
    include_internal: bool
    libclang: str | None
    regex: bool


class Decl(TypedDict):
    """One declaration found in a header."""

    name: str
    header: str
    offset: int
    kind: str
    namespace: str
    signature: str
    brief: str
    usr: str
    definition: bool


class Job(NamedTuple):
    """A translation unit to parse; see walk() for only."""

    main_file: str
    args: list[str]
    only: str | None


class TUResult(TypedDict):
    """Declarations and reached headers of one parsed translation unit."""

    file: str
    decls: list[Decl]
    seen: list[str]
    fatal: list[str]


def processCLIArgs() -> Args:
    """
    Define and process the command line interface to the runClangd.py script.
    """
    cli_description = (
        "Generate compile_commands.json and an API catalog for clangd and\n"
        "AI coding agents. All generated files are gitignored."
    )
    cli_epilog = (
        "Examples:\n"
        "  ./runClangd.py all            generate everything\n"
        "  ./runClangd.py --clean        remove everything this script made\n"
        "  ./runClangd.py --clean all    remove, then regenerate\n"
        "  ./runClangd.py cdb --opencl   OpenCL flags for every entry\n"
    )

    common = ArgumentParser(add_help=False)
    common.add_argument(
        "--compiler",
        default="clang++",
        help="clang driver used for compile_commands.json and for libclang's\n"
        "resource dir (default: clang++)",
    )
    common.add_argument(
        "--pin-system-includes",
        action="store_true",
        help="add the system include dirs reported by the compiler as\n"
        "explicit -isystem flags, for when clangd or libclang pick up\n"
        "a different standard library than the compiler does",
    )

    cdb_opts = ArgumentParser(add_help=False)
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

    cat_opts = ArgumentParser(add_help=False)
    cat_opts.add_argument(
        "-j",
        metavar="N",
        type=int,
        default=MAX_JOBS,
        help=f"number of parallel parses (max {MAX_JOBS})",
    )
    cat_opts.add_argument(
        "--out",
        default=CATALOG_DIR,
        help=f"catalog output directory (default: {CATALOG_DIR})",
    )
    cat_opts.add_argument(
        "--include-internal",
        action="store_true",
        help="list internal:: symbols inline instead of in a trailing section",
    )
    cat_opts.add_argument(
        "--libclang",
        metavar="PATH",
        help="libclang shared library, or the directory holding the clang\n"
        "python bindings (the one containing clang/cindex.py)",
    )
    cat_opts.add_argument(
        "--regex",
        action="store_true",
        help="skip libclang and use the regex fallback scanner",
    )

    parser = ArgumentParser(
        description=cli_description,
        epilog=cli_epilog,
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="remove every file this script generates (catalog at its\n"
        "default location), then run the command if one is given",
    )
    sub = parser.add_subparsers(dest="command")
    for name, parents, text in (
        ("cdb", [common, cdb_opts], "write compile_commands.json"),
        ("catalog", [common, cat_opts], "write the API catalog"),
        (
            "all",
            [common, cdb_opts, cat_opts],
            "write compile_commands.json and the API catalog",
        ),
    ):
        sub.add_parser(
            name, parents=parents, help=text, formatter_class=RawTextHelpFormatter
        )
    args = parser.parse_args(namespace=Args())
    if args.command is None and not args.clean:
        parser.print_help()
        sys.exit(1)
    if hasattr(args, "j"):
        if args.j < 1:
            stopErr("-j must be at least 1", 1)
        if args.j > MAX_JOBS:
            print(f"capping -j {args.j} at {MAX_JOBS}")
            args.j = MAX_JOBS
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
    agents = os.path.join(ROOT, ".agents")
    if os.path.isdir(agents) and not os.listdir(agents):
        print("removing .agents")
        os.rmdir(agents)


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


@functools.cache
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


############################################################
#
# API catalog
#
############################################################


def load_clang(libclang: str | None = None) -> CIndex | None:
    """Import clang.cindex and check that libclang loads, or return None."""
    if libclang and os.path.isdir(libclang):
        sys.path.insert(0, libclang)
    # The bindings live next to libclang and are often not importable until
    # the llvm-config path below is added, so pyright may not resolve them.
    try:
        from clang import cindex  # pyright: ignore[reportMissingImports]
    except ImportError:
        prefix = run_command(
            ["llvm-config", "--prefix"], capture=True, check=False
        ).stdout.strip()
        if not prefix:
            return None
        sys.path.insert(0, os.path.join(prefix, "lib", "python3", "site-packages"))
        try:
            from clang import cindex  # pyright: ignore[reportMissingImports]
        except ImportError:
            return None
    if libclang and os.path.isfile(libclang):
        cindex.Config.set_library_file(libclang)
    try:
        cindex.Index.create()
    except cindex.LibclangError:
        libdir = run_command(
            ["llvm-config", "--libdir"], capture=True, check=False
        ).stdout.strip()
        if not libdir:
            return None
        cindex.Config.loaded = False
        cindex.Config.set_library_path(libdir)
        try:
            cindex.Index.create()
        except cindex.LibclangError:
            return None
    return cindex


@dataclass
class Worker:
    """Per-process libclang state, set up by init_worker()."""

    cindex: CIndex
    index: Any
    rel: dict[str, str | None] = field(default_factory=dict)
    text: dict[str, str] = field(default_factory=dict)


_worker: Worker | None = None


def init_worker(libclang: str | None) -> None:
    """Load libclang once per worker process."""
    global _worker
    cindex = load_clang(libclang)
    if cindex is None:
        raise RuntimeError("libclang failed to load in a worker process")
    _worker = Worker(cindex, cindex.Index.create())


def worker() -> Worker:
    """This process's Worker."""
    if _worker is None:
        raise RuntimeError("init_worker() has not run in this process")
    return _worker


def rel_in_stan(name: str) -> str | None:
    """
    Repo-relative '/' path for a file libclang reports, or None if it is not
    under stan/. libclang reports -I . includes relative to the repo root.
    """
    cache = worker().rel
    if name not in cache:
        rel = os.path.relpath(os.path.join(ROOT, name), ROOT).replace(os.sep, "/")
        cache[name] = rel if rel.startswith("stan/") else None
    return cache[name]


def clip(text: str | None, width: int) -> str:
    """Collapse whitespace and shorten text to width for a one-line entry."""
    return textwrap.shorten(text or "", width, placeholder=" ...")


def file_text(name: str) -> str:
    """
    Contents of a file libclang reports, cached per worker. Decoded as
    latin-1 so each byte is one character and libclang's byte offsets index
    the text directly.
    """
    cache = worker().text
    if name not in cache:
        with open(os.path.join(ROOT, name), "rb") as f:
            cache[name] = f.read().decode("latin-1")
    return cache[name]


def scan_to_body(text: str, i: int, end: int) -> int:
    """
    Index of the first '{' or ';' in text[i:end] outside parentheses and
    brackets, or end if there is none.
    """
    depth = 0
    while i < end:
        c = text[i]
        if c in "([":
            depth += 1
        elif c in ")]":
            depth -= 1
        elif depth <= 0 and c in "{;":
            return i
        i += 1
    return end


def signature_of(cursor: Cursor) -> tuple[str, bool]:
    """
    Declaration text up to the first '{' or ';' at depth 0, and whether it
    stopped at '{'. PARSE_SKIP_FUNCTION_BODIES makes is_definition() false
    for every function and ends a definition's extent just before its '{',
    so the scan reads one byte past the extent. Source text rather than
    cursor.get_tokens() keeps macro invocations such as
    ADD_UNARY_FUNCTION(acos) instead of the macro's body.
    """
    start = cursor.extent.start
    loc = cursor.location
    if start.file is None or loc.file is None or start.file.name != loc.file.name:
        return cursor.displayname, False
    text = file_text(start.file.name)
    end = min(len(text), cursor.extent.end.offset + 1)
    i = scan_to_body(text, start.offset, end)
    has_body = text[i : i + 1] == "{"
    sig = text[start.offset : i].encode("latin-1").decode("utf-8", "replace")
    sig = re.sub(r"//[^\n]*|/\*.*?\*/", " ", sig, flags=re.DOTALL)
    return clip(sig, 400) or cursor.displayname, has_body


def first_sentence(text: str | None) -> str:
    """First sentence of a doc comment, with comment markers removed."""
    if not text:
        return ""
    text = re.sub(r"^\s*/\*[*!]?|\*/\s*$", "", text)
    lines: list[str] = []
    for line in text.splitlines():
        line = re.sub(r"^\s*(\*|///|//!)\s?", "", line).strip()
        if not line:
            if lines:
                break
            continue
        if line.startswith(("@", "\\")) and not re.match(r"[@\\]brief\b", line):
            break
        lines.append(line)
    text = re.sub(r"^[@\\]brief\s*", "", " ".join(lines))
    match = re.match(r"(.+?[.!?])(\s|$)", text)
    return match.group(1) if match else text


def walk(
    cursor: Cursor, namespaces: list[str], only: str | None, out: list[Decl]
) -> None:
    """
    Collect declarations under a namespace cursor, recursing into namespaces.
    With only set, keep declarations whose header starts with it.
    """
    for child in cursor.get_children():
        loc = child.location
        if loc.file is None:
            continue
        rel = rel_in_stan(loc.file.name)
        if rel is None or (only and not rel.startswith(only)):
            continue
        kind_name: str = child.kind.name
        if kind_name == "NAMESPACE":
            walk(child, [*namespaces, child.spelling], only, out)
            continue
        kind = KEEP_KINDS.get(kind_name)
        if kind is None or not child.spelling:
            continue
        if kind_name == "VAR_DECL":
            parent = child.semantic_parent
            if parent is not None and parent.kind.name != "NAMESPACE":
                continue
        signature, has_body = signature_of(child)
        out.append(
            Decl(
                name=child.spelling,
                header=rel,
                offset=loc.offset,
                kind=kind,
                namespace="::".join(namespaces),
                signature=signature,
                brief=clip(
                    child.brief_comment or first_sentence(child.raw_comment), 240
                ),
                usr=child.get_usr(),
                definition=has_body
                or child.is_definition()
                or kind in ("alias", "variable"),
            )
        )


def parse_tu(job: Job) -> TUResult:
    """Parse one translation unit and return its declarations."""
    w = worker()
    cindex = w.cindex
    opts = (
        cindex.TranslationUnit.PARSE_SKIP_FUNCTION_BODIES
        | cindex.TranslationUnit.PARSE_INCOMPLETE
    )
    try:
        tu = w.index.parse(
            os.path.join(ROOT, job.main_file), args=job.args, options=opts
        )
    except cindex.TranslationUnitLoadError as e:
        return TUResult(file=job.main_file, decls=[], seen=[], fatal=[str(e)])
    fatal = [str(d) for d in tu.diagnostics if d.severity >= cindex.Diagnostic.Fatal]
    seen = {job.main_file}
    for inc in tu.get_includes():
        rel = rel_in_stan(inc.include.name)
        if rel is not None:
            seen.add(rel)
    decls: list[Decl] = []
    for child in tu.cursor.get_children():
        if child.kind == cindex.CursorKind.NAMESPACE and child.spelling == "stan":
            walk(child, ["stan"], job.only, decls)
    return TUResult(file=job.main_file, decls=decls, seen=sorted(seen), fatal=fatal)


def module_of(header: str) -> str:
    """Catalog file stem for a header, e.g. prim-fun for stan/math/prim/fun/x.hpp."""
    parts = header.split("/")
    if len(parts) >= 5:
        return parts[2] + "-" + parts[3]
    if len(parts) == 4:
        return parts[2]
    return "math"


def display_name(decl: Decl) -> str:
    """Name as grepped in the catalog: qualified only outside stan::math."""
    ns = decl["namespace"].split("::")
    ns = [n for n in ns if n and n != "internal"]
    if ns[:2] == ["stan", "math"]:
        ns = ns[2:]
    elif ns[:1] == ["stan"]:
        ns = ns[1:]
    return "::".join([*ns, decl["name"]])


def catalog_libclang(args: Args) -> list[Decl]:
    """Collect declarations for every header with libclang."""
    # Without the compiler's resource dir libclang may miss its builtin
    # headers (stddef.h) and silently turn unknown types into int.
    extra = [
        "-resource-dir",
        run_command(
            [args.compiler, "-print-resource-dir"], capture=True
        ).stdout.strip(),
        "-Wno-unknown-warning-option",
        "-w",
    ]
    if args.pin_system_includes:
        extra += system_includes(args.compiler)

    def tu_args(opencl: bool) -> list[str]:
        return ["-x", "c++", *build_flags(opencl)[0], *extra]

    def run_pass(
        pool: Executor, label: str, jobs: list[Job]
    ) -> tuple[list[TUResult], list[tuple[str, str]]]:
        t0 = time.time()
        chunksize = max(1, len(jobs) // (4 * args.j))
        results = list(pool.map(parse_tu, jobs, chunksize=chunksize))
        print(f"{label} pass: {len(jobs)} TUs ({time.time() - t0:.1f}s)")
        return results, [(r["file"], r["fatal"][0]) for r in results if r["fatal"]]

    # rev and prim come from mix.hpp, so the OpenCL TU only keeps opencl/.
    umbrella = [
        Job("stan/math/mix.hpp", tu_args(False), None),
        Job("stan/math/opencl/rev.hpp", tu_args(True), "stan/math/opencl/"),
    ]
    headers = find_files("stan", ".hpp")
    with ProcessPoolExecutor(
        max_workers=args.j, initializer=init_worker, initargs=(args.libclang,)
    ) as pool:
        results, failed = run_pass(pool, "umbrella", umbrella)
        if failed:
            file, msg = failed[0]
            stopErr(
                f"fatal errors parsing {file}:\n  {msg}\nif builtin headers are"
                " missing, pass a --compiler matching libclang's version",
                1,
            )
        seen: set[str] = set()
        for res in results:
            seen.update(res["seen"])
        stragglers = [h for h in headers if h not in seen]
        print(f"umbrella pass reached {len(seen)} of {len(headers)} headers")
        jobs = [Job(h, tu_args(is_opencl_path(h)), h) for h in stragglers]
        more, failed = run_pass(pool, "straggler", jobs)
        results += more
    if failed:
        print(
            f"{len(failed)} headers had fatal errors on their own; partial results kept:"
        )
        for file, msg in failed[:20]:
            print(f"  {file}: {msg}")
    decls: dict[tuple[str, int], Decl] = {}
    for res in results:
        for d in res["decls"]:
            decls.setdefault((d["header"], d["offset"]), d)
    return list(decls.values())


REGEX_DECL = re.compile(
    r"^[ \t]*(?:template\s*<(?P<tparams>[^;{]*?)>\s*)?"
    r"(?:(?:inline|static|constexpr|STAN_COLD_PATH|explicit)\s+)*"
    r"(?:(?P<cls>struct|class)\s+(?P<cname>\w+)"
    r"|using\s+(?P<aname>\w+)\s*="
    r"|[\w:<>,\s\*&]+?\b(?P<fname>\w+)\s*\()",
    re.MULTILINE,
)


def scopes(text: str) -> tuple[str, list[int], list[tuple[bool, bool]]]:
    """
    Blank out comments and string literals (keeping offsets), then return
    that code, the sorted offsets where a brace scan changes scope, and the
    (at_namespace_scope, in_internal) state starting at each offset.
    """
    code = re.sub(
        r"//[^\n]*|/\*.*?\*/|R\"\((.*?)\)\"|\"(\\.|[^\"\\])*\"",
        lambda m: re.sub(r"[^\n]", " ", m.group(0)),
        text,
        flags=re.DOTALL,
    )
    stack: list[str] = []
    offsets = [0]
    states = [(True, False)]
    last = 0
    for i, c in enumerate(code):
        if c == ";":
            last = i + 1
            continue
        if c == "{":
            match = re.search(r"namespace\s*(\w*)\s*$", code[last:i])
            if match is None:
                stack.append("other")
            else:
                stack.append("internal" if match.group(1) == "internal" else "ns")
        elif c == "}" and stack:
            stack.pop()
        else:
            continue
        last = i + 1
        offsets.append(i + 1)
        states.append(("other" not in stack, "internal" in stack))
    return code, offsets, states


def catalog_regex() -> list[Decl]:
    """Collect declarations with a regex scan, when libclang is unavailable."""
    decls: list[Decl] = []
    doc = re.compile(r"/\*\*(.*?)\*/\s*$", re.DOTALL)
    for rel in find_files("stan", ".hpp"):
        with open(os.path.join(ROOT, rel), errors="replace") as f:
            text = f.read()
        code, offsets, states = scopes(text)
        for m in REGEX_DECL.finditer(code):
            name = m.group("cname") or m.group("aname") or m.group("fname")
            if not name or name in ("if", "for", "while", "switch", "return"):
                continue
            at_ns, in_internal = states[bisect.bisect_right(offsets, m.start()) - 1]
            if not at_ns:
                continue
            if m.group("cname"):
                kind = "class"
            elif m.group("aname"):
                kind = "alias"
            else:
                kind = "function"
            dm = doc.search(text[max(0, m.start() - 4000) : m.start()])
            end = scan_to_body(code, m.end(), len(code))
            ns = ["stan", "math"] + (["internal"] if in_internal else [])
            brief = first_sentence("/**" + dm.group(1) + "*/") if dm else ""
            decls.append(
                Decl(
                    name=name,
                    header=rel,
                    offset=m.start(),
                    kind=kind,
                    namespace="::".join(ns),
                    signature=clip(code[m.start() : end], 400),
                    brief=clip(brief + " [regex]", 240),
                    usr="",
                    definition=True,
                )
            )
    return decls


def write_catalog(args: Args) -> None:
    """Write the API catalog to args.out."""
    out_rel = os.path.relpath(os.path.abspath(args.out), ROOT)
    require_ignored(os.path.join(out_rel, "index.md"))
    start = time.time()
    cindex = None if args.regex else load_clang(args.libclang)
    if cindex is None:
        if not args.regex:
            print("libclang not found (try --libclang PATH); using the regex scanner")
        decls = catalog_regex()
        source = "regex scan"
    else:
        decls = catalog_libclang(args)
        source = "libclang"

    # Drop forward declarations of things defined somewhere. The USR ties a
    # declaration to its definition; the regex scan has none.
    def ident(d: Decl) -> str | tuple[str, str]:
        return d["usr"] or (d["name"], d["kind"])

    defined = {ident(d) for d in decls if d["definition"]}
    decls = [d for d in decls if d["definition"] or ident(d) not in defined]
    groups: defaultdict[tuple[str, str, str, str], list[Decl]] = defaultdict(list)
    for d in sorted(decls, key=lambda d: (d["header"], d["offset"])):
        internal = "internal" in d["namespace"].split("::")
        section = "internal" if internal and not args.include_internal else "public"
        key = (module_of(d["header"]), section, display_name(d), d["header"])
        groups[key].append(d)
    lines: defaultdict[str, dict[str, list[str]]] = defaultdict(
        lambda: {"public": [], "internal": []}
    )
    index: defaultdict[str, set[str]] = defaultdict(set)
    for (module, section, name, header), ds in groups.items():
        chosen = next((d for d in ds if d["brief"]), ds[0])
        sig = chosen["signature"]
        if len(ds) > 1:
            sig += f" (x{len(ds)})"
        briefs = list(dict.fromkeys(d["brief"] for d in ds if d["brief"]))
        brief = clip(" / ".join(briefs[:3]), 240)
        lines[module][section].append(f"{name} | {header} | {sig} | {brief}")
        if section == "public":
            index[name].add(module.replace("-", "/"))
    sha = run_command(
        ["git", "rev-parse", "--short", "HEAD"], capture=True
    ).stdout.strip()
    today = datetime.datetime.now().astimezone().date().isoformat()
    stamp = (
        f"Generated by ./runClangd.py catalog ({source}) at {today} from {sha}."
        " Do not edit."
    )
    out = os.path.join(ROOT, out_rel)
    tmp = out + ".tmp"
    if os.path.isdir(tmp):
        shutil.rmtree(tmp)
    os.makedirs(tmp)
    for module in sorted(lines):
        body = [
            f"# Stan Math API catalog: {module.replace('-', '/')}",
            stamp,
            "Format: name | header | signature | brief",
            "",
        ]
        body += sorted(lines[module]["public"])
        if lines[module]["internal"]:
            body += ["", "## internal", "", *sorted(lines[module]["internal"])]
        with open(os.path.join(tmp, module + ".md"), "w") as f:
            f.write("\n".join(body) + "\n")
    with open(os.path.join(tmp, "index.md"), "w") as f:
        f.write(f"# Stan Math API catalog index\n{stamp}\n")
        f.write(
            "Format: name: modules defining it (catalog file = module with / -> -)\n\n"
        )
        for name in sorted(index):
            f.write(f"{name}: {' '.join(sorted(index[name]))}\n")
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.replace(tmp, out)
    entries = sum(len(v["public"]) + len(v["internal"]) for v in lines.values())
    print(
        f"wrote {out_rel}: {len(lines)} files, {entries} entries,"
        f" {len(index)} names in index.md ({time.time() - start:.1f}s)"
    )


def main() -> None:
    args = processCLIArgs()
    os.chdir(ROOT)
    if args.clean:
        clean()
    if args.command in ("cdb", "all"):
        write_cdb(args)
    if args.command in ("catalog", "all"):
        write_catalog(args)


if __name__ == "__main__":
    main()

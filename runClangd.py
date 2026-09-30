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

import bisect
import datetime
import functools
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import textwrap
import time
from argparse import ArgumentParser, RawTextHelpFormatter
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

ROOT = os.path.dirname(os.path.abspath(__file__))
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


def processCLIArgs():
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
        help="number of parallel parses (max %d)" % MAX_JOBS,
    )
    cat_opts.add_argument(
        "--out",
        default=CATALOG_DIR,
        help="catalog output directory (default: %s)" % CATALOG_DIR,
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
    args = parser.parse_args()
    if args.command is None and not args.clean:
        parser.print_help()
        sys.exit(1)
    if hasattr(args, "j"):
        if args.j < 1:
            stopErr("-j must be at least 1", 1)
        if args.j > MAX_JOBS:
            print("capping -j %d at %d" % (args.j, MAX_JOBS))
            args.j = MAX_JOBS
    return args


def stopErr(msg, returncode):
    """Report an error message to stderr and exit with a given code."""
    sys.stderr.write("%s\n" % msg)
    sys.stderr.write("exit now (%s)\n" % time.strftime("%x %X %Z"))
    sys.exit(returncode)


def run(command):
    """Run a command in the repo root and return its stdout; exit on failure."""
    try:
        proc = subprocess.run(
            command,
            cwd=ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
    except OSError as e:
        stopErr("cannot run %s: %s" % (command[0], e), 1)
    if proc.returncode != 0:
        stopErr(
            "command failed: %s\n%s" % (" ".join(command), proc.stderr),
            proc.returncode,
        )
    return proc.stdout


def git_ok(*args):
    """True if the git command succeeds."""
    proc = subprocess.run(
        ["git"] + list(args),
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return proc.returncode == 0


def require_ignored(path):
    """Abort unless git ignores path, so generated files are never committed."""
    if not git_ok("check-ignore", "-q", path):
        stopErr("%s is not gitignored; refusing to generate it" % path, 1)


def clean():
    """Remove every path in GENERATED_PATHS that exists and git does not track."""
    for rel in GENERATED_PATHS:
        path = os.path.join(ROOT, rel)
        if not os.path.lexists(path):
            continue
        if git_ok("ls-files", "--error-unmatch", rel):
            print("refusing to remove %s: it is tracked by git" % rel)
            continue
        print("removing %s" % rel)
        if os.path.isdir(path) and not os.path.islink(path):
            shutil.rmtree(path)
        else:
            os.remove(path)
    agents = os.path.join(ROOT, ".agents")
    if os.path.isdir(agents) and not os.listdir(agents):
        print("removing .agents")
        os.rmdir(agents)


def make_vars(names, opencl=False):
    """Read make variables through the print-% rule, so make/local is honored."""
    command = ["make", "-s"]
    if opencl:
        command.append("STAN_OPENCL=true")
    command += ["print-" + name for name in names]
    values = {}
    for line in run(command).splitlines():
        match = re.match(r"^(\w+) = ?(.*)$", line)
        if match and match.group(1) in names:
            values[match.group(1)] = shlex.split(match.group(2))
    missing = [name for name in names if name not in values]
    if missing:
        stopErr("make did not print: %s" % ", ".join(missing), 1)
    return values


@functools.lru_cache(maxsize=None)
def build_flags(opencl=False):
    """Compiler flags for Stan headers and for gtest files."""
    v = make_vars(
        ["CXXFLAGS", "CPPFLAGS", "INC_GTEST", "CXXFLAGS_GTEST", "CPPFLAGS_GTEST"],
        opencl=opencl,
    )
    base = v["CXXFLAGS"] + v["CPPFLAGS"]
    gtest = v["INC_GTEST"] + v["CXXFLAGS_GTEST"] + v["CPPFLAGS_GTEST"]
    return base, gtest


@functools.lru_cache(maxsize=None)
def system_includes(compiler):
    """-isystem flags for the dirs the compiler searches for <...> includes."""
    proc = subprocess.run(
        [compiler, "-E", "-v", "-x", "c++", os.devnull],
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    flags = []
    inside = False
    for line in proc.stderr.splitlines():
        if line.startswith("#include <...> search starts here"):
            inside = True
        elif line.startswith("End of search list"):
            break
        elif inside:
            flags += ["-isystem", os.path.normpath(line.strip())]
    return flags


def find_files(top, suffix):
    """Sorted repo-relative '/' paths under top ending in suffix."""
    found = []
    for dirpath, _, filenames in os.walk(os.path.join(ROOT, top)):
        rel_dir = os.path.relpath(dirpath, ROOT).replace(os.sep, "/")
        found += [rel_dir + "/" + name for name in filenames if name.endswith(suffix)]
    return sorted(found)


def is_opencl_path(rel):
    """True for files that only compile with STAN_OPENCL defined."""
    return rel.startswith(("stan/math/opencl/", "test/unit/math/opencl/"))


def write_cdb(args):
    """Write compile_commands.json for every header and unit test."""
    require_ignored(CDB_FILE)
    start = time.time()
    pinned = system_includes(args.compiler) if args.pin_system_includes else []

    def entry(rel, head, tail):
        base, gtest = build_flags(args.opencl or is_opencl_path(rel))
        flags = base + (gtest if tail[0] == "-c" else [])
        return {
            "directory": ROOT,
            "file": os.path.join(ROOT, rel),
            "arguments": [args.compiler] + head + flags + pinned + tail,
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
        "wrote %s: %d headers, %d tests (%.1fs)"
        % (CDB_FILE, len(headers), len(tests), time.time() - start)
    )
    print("rerun after changing make/local or adding files.")


############################################################
#
# API catalog
#
############################################################


def llvm_config(flag):
    """Output of llvm-config FLAG, or '' if llvm-config is unavailable."""
    try:
        return subprocess.run(
            ["llvm-config", flag],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            universal_newlines=True,
        ).stdout.strip()
    except OSError:
        return ""


def load_clang(libclang=None):
    """Import clang.cindex and check that libclang loads, or return None."""
    if libclang and os.path.isdir(libclang):
        sys.path.insert(0, libclang)
    try:
        import clang.cindex as cindex
    except ImportError:
        prefix = llvm_config("--prefix")
        if not prefix:
            return None
        sys.path.insert(0, os.path.join(prefix, "lib", "python3", "site-packages"))
        try:
            import clang.cindex as cindex
        except ImportError:
            return None
    if libclang and os.path.isfile(libclang):
        cindex.Config.set_library_file(libclang)
    try:
        cindex.Index.create()
    except Exception:
        libdir = llvm_config("--libdir")
        if not libdir:
            return None
        cindex.Config.loaded = False
        cindex.Config.set_library_path(libdir)
        try:
            cindex.Index.create()
        except Exception:
            return None
    return cindex


_worker = {}


def init_worker(libclang):
    """Load libclang once per worker process."""
    cindex = load_clang(libclang)
    _worker["cindex"] = cindex
    _worker["index"] = cindex.Index.create()
    _worker["rel"] = {}
    _worker["text"] = {}


def rel_in_stan(name):
    """
    Repo-relative '/' path for a file libclang reports, or None if it is not
    under stan/. libclang reports -I . includes relative to the repo root.
    """
    cache = _worker["rel"]
    if name not in cache:
        rel = os.path.relpath(os.path.join(ROOT, name), ROOT).replace(os.sep, "/")
        cache[name] = rel if rel.startswith("stan/") else None
    return cache[name]


def clip(text, width):
    """Collapse whitespace and shorten text to width for a one-line entry."""
    return textwrap.shorten(text or "", width, placeholder=" ...")


def file_bytes(name):
    """Contents of a file libclang reports, cached per worker."""
    cache = _worker["text"]
    if name not in cache:
        with open(os.path.join(ROOT, name), "rb") as f:
            cache[name] = f.read()
    return cache[name]


def signature_of(cursor):
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
    text = file_bytes(start.file.name)
    depth = 0
    i = start.offset
    end = min(len(text), cursor.extent.end.offset + 1)
    while i < end:
        c = text[i : i + 1]
        if c in (b"(", b"["):
            depth += 1
        elif c in (b")", b"]"):
            depth -= 1
        elif depth == 0 and c in (b"{", b";"):
            break
        i += 1
    has_body = text[i : i + 1] == b"{"
    sig = text[start.offset : i].decode("utf-8", "replace")
    sig = re.sub(r"//[^\n]*|/\*.*?\*/", " ", sig, flags=re.S)
    return clip(sig, 400) or cursor.displayname, has_body


def first_sentence(text):
    """First sentence of a doc comment, with comment markers removed."""
    if not text:
        return ""
    text = re.sub(r"^\s*/\*[*!]?|\*/\s*$", "", text)
    lines = []
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


def walk(cursor, namespaces, only, out):
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
        kind_name = child.kind.name
        if kind_name == "NAMESPACE":
            walk(child, namespaces + [child.spelling], only, out)
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
            {
                "name": child.spelling,
                "header": rel,
                "offset": loc.offset,
                "kind": kind,
                "namespace": "::".join(namespaces),
                "signature": signature,
                "brief": clip(
                    child.brief_comment or first_sentence(child.raw_comment), 240
                ),
                "usr": child.get_usr(),
                "definition": has_body
                or child.is_definition()
                or kind in ("alias", "variable"),
            }
        )


def parse_tu(job):
    """
    Parse one translation unit and return its declarations.

    job is (main_file, args, only); see walk() for only.
    """
    main_file, args, only = job
    cindex = _worker["cindex"]
    opts = (
        cindex.TranslationUnit.PARSE_SKIP_FUNCTION_BODIES
        | cindex.TranslationUnit.PARSE_INCOMPLETE
    )
    try:
        tu = _worker["index"].parse(
            os.path.join(ROOT, main_file), args=args, options=opts
        )
    except cindex.TranslationUnitLoadError as e:
        return {"file": main_file, "decls": [], "seen": [], "fatal": [str(e)]}
    fatal = [str(d) for d in tu.diagnostics if d.severity >= cindex.Diagnostic.Fatal]
    seen = {main_file}
    for inc in tu.get_includes():
        rel = rel_in_stan(inc.include.name)
        if rel is not None:
            seen.add(rel)
    decls = []
    for child in tu.cursor.get_children():
        if child.kind == cindex.CursorKind.NAMESPACE and child.spelling == "stan":
            walk(child, ["stan"], only, decls)
    return {"file": main_file, "decls": decls, "seen": sorted(seen), "fatal": fatal}


def module_of(header):
    """Catalog file stem for a header, e.g. prim-fun for stan/math/prim/fun/x.hpp."""
    parts = header.split("/")
    if len(parts) >= 5:
        return parts[2] + "-" + parts[3]
    if len(parts) == 4:
        return parts[2]
    return "math"


def display_name(decl):
    """Name as grepped in the catalog: qualified only outside stan::math."""
    ns = decl["namespace"].split("::")
    ns = [n for n in ns if n and n != "internal"]
    if ns[:2] == ["stan", "math"]:
        ns = ns[2:]
    elif ns[:1] == ["stan"]:
        ns = ns[1:]
    return "::".join(ns + [decl["name"]])


def catalog_libclang(args):
    """Collect declarations for every header with libclang."""
    # Without the compiler's resource dir libclang may miss its builtin
    # headers (stddef.h) and silently turn unknown types into int.
    extra = [
        "-resource-dir",
        run([args.compiler, "-print-resource-dir"]).strip(),
        "-Wno-unknown-warning-option",
        "-w",
    ]
    if args.pin_system_includes:
        extra += system_includes(args.compiler)

    def tu_args(opencl):
        return ["-x", "c++"] + build_flags(opencl)[0] + extra

    def run_pass(pool, label, jobs):
        t0 = time.time()
        chunksize = max(1, len(jobs) // (4 * args.j))
        results = list(pool.map(parse_tu, jobs, chunksize=chunksize))
        print("%s pass: %d TUs (%.1fs)" % (label, len(jobs), time.time() - t0))
        return results, [(r["file"], r["fatal"][0]) for r in results if r["fatal"]]

    # rev and prim come from mix.hpp, so the OpenCL TU only keeps opencl/.
    umbrella = [
        ("stan/math/mix.hpp", tu_args(False), None),
        ("stan/math/opencl/rev.hpp", tu_args(True), "stan/math/opencl/"),
    ]
    headers = find_files("stan", ".hpp")
    with ProcessPoolExecutor(
        max_workers=args.j, initializer=init_worker, initargs=(args.libclang,)
    ) as pool:
        results, failed = run_pass(pool, "umbrella", umbrella)
        if failed:
            stopErr(
                "fatal errors parsing %s:\n  %s\nif builtin headers are missing,"
                " pass a --compiler matching libclang's version" % failed[0],
                1,
            )
        seen = set()
        for res in results:
            seen.update(res["seen"])
        stragglers = [h for h in headers if h not in seen]
        print("umbrella pass reached %d of %d headers" % (len(seen), len(headers)))
        jobs = [(h, tu_args(is_opencl_path(h)), h) for h in stragglers]
        more, failed = run_pass(pool, "straggler", jobs)
        results += more
    if failed:
        print("%d headers had fatal errors on their own; partial results kept:" % len(failed))
        for f, msg in failed[:20]:
            print("  %s: %s" % (f, msg))
    decls = {}
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
    re.M,
)


def scopes(text):
    """
    Blank out comments and string literals (keeping offsets), then return
    that code, the sorted offsets where a brace scan changes scope, and the
    (at_namespace_scope, in_internal) state starting at each offset.
    """
    code = re.sub(
        r"//[^\n]*|/\*.*?\*/|R\"\((.*?)\)\"|\"(\\.|[^\"\\])*\"",
        lambda m: re.sub(r"[^\n]", " ", m.group(0)),
        text,
        flags=re.S,
    )
    stack = []
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


def catalog_regex():
    """Collect declarations with a regex scan, when libclang is unavailable."""
    decls = []
    doc = re.compile(r"/\*\*(.*?)\*/\s*$", re.S)
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
            kind = "class" if m.group("cname") else "alias" if m.group("aname") else "function"
            dm = doc.search(text[max(0, m.start() - 4000) : m.start()])
            end = m.end()
            depth = 0
            while end < len(code):
                c = code[end]
                if c in "([":
                    depth += 1
                elif c in ")]":
                    depth -= 1
                elif depth <= 0 and c in "{;":
                    break
                end += 1
            ns = ["stan", "math"] + (["internal"] if in_internal else [])
            brief = first_sentence("/**" + dm.group(1) + "*/") if dm else ""
            decls.append(
                {
                    "name": name,
                    "header": rel,
                    "offset": m.start(),
                    "kind": kind,
                    "namespace": "::".join(ns),
                    "signature": clip(code[m.start() : end], 400),
                    "brief": clip(brief + " [regex]", 240),
                    "usr": "",
                    "definition": True,
                }
            )
    return decls


def write_catalog(args):
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
    def ident(d):
        return d["usr"] or (d["name"], d["kind"])

    defined = {ident(d) for d in decls if d["definition"]}
    decls = [d for d in decls if d["definition"] or ident(d) not in defined]
    groups = defaultdict(list)
    for d in sorted(decls, key=lambda d: (d["header"], d["offset"])):
        internal = "internal" in d["namespace"].split("::")
        section = "internal" if internal and not args.include_internal else "public"
        key = (module_of(d["header"]), section, display_name(d), d["header"])
        groups[key].append(d)
    lines = defaultdict(lambda: {"public": [], "internal": []})
    index = defaultdict(set)
    for (module, section, name, header), ds in groups.items():
        chosen = next((d for d in ds if d["brief"]), ds[0])
        sig = chosen["signature"]
        if len(ds) > 1:
            sig += " (x%d)" % len(ds)
        briefs = list(dict.fromkeys(d["brief"] for d in ds if d["brief"]))
        brief = clip(" / ".join(briefs[:3]), 240)
        lines[module][section].append("%s | %s | %s | %s" % (name, header, sig, brief))
        if section == "public":
            index[name].add(module.replace("-", "/"))
    sha = run(["git", "rev-parse", "--short", "HEAD"]).strip()
    stamp = "Generated by ./runClangd.py catalog (%s) at %s from %s. Do not edit." % (
        source,
        datetime.date.today().isoformat(),
        sha,
    )
    out = os.path.join(ROOT, out_rel)
    tmp = out + ".tmp"
    if os.path.isdir(tmp):
        shutil.rmtree(tmp)
    os.makedirs(tmp)
    for module in sorted(lines):
        body = [
            "# Stan Math API catalog: %s" % module.replace("-", "/"),
            stamp,
            "Format: name | header | signature | brief",
            "",
        ]
        body += sorted(lines[module]["public"])
        if lines[module]["internal"]:
            body += ["", "## internal", ""] + sorted(lines[module]["internal"])
        with open(os.path.join(tmp, module + ".md"), "w") as f:
            f.write("\n".join(body) + "\n")
    with open(os.path.join(tmp, "index.md"), "w") as f:
        f.write("# Stan Math API catalog index\n%s\n" % stamp)
        f.write("Format: name: modules defining it (catalog file = module with / -> -)\n\n")
        for name in sorted(index):
            f.write("%s: %s\n" % (name, " ".join(sorted(index[name]))))
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.replace(tmp, out)
    entries = sum(len(v["public"]) + len(v["internal"]) for v in lines.values())
    print(
        "wrote %s: %d files, %d entries, %d names in index.md (%.1fs)"
        % (out_rel, len(lines), entries, len(index), time.time() - start)
    )


def main():
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

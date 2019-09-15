#!/usr/bin/python

"""
Replacement for test-math-dependencies target in Makefile.

"""

from __future__ import print_function
import os
import sys
import re
import glob

winsfx = ".exe"
testsfx = "_test.cpp"


def files_in_folder(folder):
    """Returns a list of files in the folder and all
    its subfolders recursively. The folder can be
    written with wildcards as with the Unix find command.
    """
    files = []
    for f in glob.glob(folder):
        if os.path.isdir(f):
            files.extend(files_in_folder(f + os.sep + "**"))
        else:
            files.append(f)
    return files


def grep_patterns(type, folder, patterns_and_messages, exclude_filters=[]):
    """Checks the files in the provided folder for matches
    with any of the patterns. It returns an array of
    messages and the provided type with
    the line number. This check ignores comments.
    @param type: type or group of the check, listed with the error
    @param folder: folder in which to check for the pattern
    @param patterns_and_messages: a list of patterns and messages that 
        are printed if the pattern is matched
    @param exclude_filter a list of files or folder that are excluded from
        the check
    """
    errors = []
    folder.replace("/", os.sep)
    exclude_files = []
    for excl in exclude_filters:
        exclude_files.extend(files_in_folder(excl))
    files = files_in_folder(folder + os.sep + "**")
    files = [x for x in files if x not in exclude_files]
    for filepath in files:
        if os.path.isfile(filepath):
            line_num = 0
            multi_line_comment = False
            old_state_multi_line_comment = False
            with open(filepath, "r") as f:
                for line in f:
                    line_num += 1
                    # exclude multi line comments
                    if multi_line_comment:
                        if re.search("\*/", line):
                            multi_line_comment = False
                    else:
                        if re.search("/\*", line):
                            multi_line_comment = True
                    # parse the first line in a multi line comment for rare and weird case of
                    # "pattern /*""
                    if not multi_line_comment or (
                        multi_line_comment and not old_state_multi_line_comment
                    ):
                        for p in patterns_and_messages:
                            # cover the edge cases where matched patterns
                            # are behind "//", "/*" or before "*/"
                            if (
                                not re.search(".*" + p["pattern"] + ".*\*/.*", line)
                                and not re.search(".*/\*.*" + p["pattern"], line)
                                and not re.search(".*//.*" + p["pattern"], line)
                                and re.search(p["pattern"], line)
                            ):
                                errors.append(
                                    filepath
                                    + " at line "
                                    + str(line_num)
                                    + ":\n\t"
                                    + "["
                                    + type
                                    + "] "
                                    + p["message"]
                                )
                    old_state_multi_line_comment = multi_line_comment
    return errors


def check_non_test_files_in_test():
    all_cpp = files_in_folder("test/unit/math/")
    # if the file is a .cpp file that doesnt end with _test.cpp
    errors = [
        x + ":\n\t A .cpp file without the " + testsfx + " suffix in test/unit/math/"
        for x in all_cpp
        if os.path.splitext(x)[1] == ".cpp" and x[-len(testsfx) :] != testsfx
    ]
    return errors


def main():
    errors = []
    # Check for files inside stan/math/prim that contain stan/math/rev or stan/math/fwd
    prim_checks = [
        {
            "pattern": "<stan/math/rev/",
            "message": "File includes a stan/math/rev header file.",
        },
        {"pattern": "stan::math::var", "message": "File uses stan::math::var."},
        {
            "pattern": "<stan/math/fwd/",
            "message": "File includes a stan/math/fwd header file.",
        },
        {"pattern": "stan::math::fvar", "message": "File uses stan::math::fvar."},
        {
            "pattern": "<stan/math/mix/",
            "message": "File includes a stan/math/mix header file.",
        },
    ]
    errors.extend(grep_patterns("prim", "stan/math/prim", prim_checks))

    # Check for files inside stan/math/rev that contain stan/math/fwd or stan/math/mix
    rev_checks = [
        {
            "pattern": "<stan/math/fwd/",
            "message": "File includes a stan/math/fwd header file.",
        },
        {"pattern": "stan::math::fvar", "message": "File uses stan::math::fvar."},
        {
            "pattern": "<stan/math/mix/",
            "message": "File includes a stan/math/mix header file.",
        },
    ]
    errors.extend(grep_patterns("rev", "stan/math/rev", rev_checks))

    # Check for files inside stan/math/*/scal that contain stan/math/*/arr or stan/math/*/mat
    scal_checks = [
        {
            "pattern": "<stan/math/.*/arr/",
            "message": "File includes an array header file.",
        },
        {"pattern": "<vector>", "message": "File includes an std::vector header."},
        {"pattern": "std::vector", "message": "File uses std::vector."},
        {
            "pattern": "<stan/math/.*/mat/",
            "message": "File includes a matrix header file.",
        },
        {"pattern": "<Eigen", "message": "File includes an Eigen header."},
        {"pattern": "Eigen::", "message": "File uses Eigen."},
    ]
    errors.extend(grep_patterns("scal", "stan/math/*/scal", scal_checks))

    # Check for files inside stan/math/*/arr that contain stan/math/*/mat or Eigen
    arr_checks = [
        {
            "pattern": "<stan/math/.*/mat/",
            "message": "File includes an matrix header file.",
        },
        {"pattern": "<Eigen", "message": "File includes an Eigen header."},
        {"pattern": "Eigen::", "message": "File uses Eigen."},
    ]
    errors.extend(grep_patterns("arr", "stan/math/*/arr", arr_checks))

    # Check to make sure we use C++14 constructs in stan/math
    cpp14_checks = [
        {
            "pattern": "boost::is_unsigned<",
            "message": "File uses boost::is_unsigned instead of std::is_unsigned.",
        },
        {
            "pattern": "<boost/type_traits/is_unsigned>",
            "message": "File includes <boost/type_traits/is_unsigned.hpp> instead of <type_traits>.",
        },
        {
            "pattern": "boost::is_arithmetic<",
            "message": "File uses boost::is_arithmetic instead of std::is_arithmetic.",
        },
        {
            "pattern": "<boost/type_traits/is_arithmetic.hpp>",
            "message": "File includes <boost/type_traits/is_unsigned.hpp> instead of <type_traits>.",
        },
        {
            "pattern": "boost::is_convertible<",
            "message": "File uses boost::is_convertible instead of std::is_convertible.",
        },
        {
            "pattern": "<boost/type_traits/is_convertible.hpp>",
            "message": "File includes <boost/type_traits/is_convertible.hpp> instead of <type_traits>.",
        },
        {
            "pattern": "boost::is_same<",
            "message": "File uses boost::is_same instead of std::is_same.",
        },
        {
            "pattern": "<boost/type_traits/is_same.hpp>",
            "message": "File includes <boost/type_traits/is_same.hpp> instead of <type_traits>.",
        },
        {
            "pattern": "boost::enable_if_c<",
            "message": "File uses boost::enable_if_c instead of std::enable_if.",
        },
        {
            "pattern": "boost::enable_if<",
            "message": "File uses boost::enable_if instead of std::enable_if.",
        },
        {
            "pattern": "boost::disable_if<",
            "message": "File uses boost::disable_if instead of std::enable_if.",
        },
        {
            "pattern": "<boost/utility/enable_if.hpp>",
            "message": "Replace \<boost/utility/enable_if.hpp\> with \<type_traits\>.",
        },
    ]
    errors.extend(grep_patterns("C++14", "stan/math", cpp14_checks))

    # Check for includes of stan/math/*/meta/*.hpp inside stan/math, excluding meta.hpp files and the /meta subfolder
    meta_checks = [
        {
            "pattern": "<stan/math/.*/meta/.*hpp",
            "message": "File includes */meta/*.hpp header file. Should include meta.hpp",
        }
    ]
    meta_exclude = ["stan/math/*/*/meta", "stan/math/*/meta.hpp"]
    errors.extend(grep_patterns("meta", "stan/math", meta_checks, meta_exclude))

    #  Check that we do not use non-reentrant safe functions from std
    thread_safe_checks = [
        {
            "pattern": "std::lgamma",
            "message": "Replace std::lgamma with reentrant safe version boost::math::lgamma or lgamma_r on supporting platforms.",
        }
    ]
    errors.extend(
        grep_patterns("thread-reentrant-safe", "stan/math", thread_safe_checks)
    )

    errors.extend(check_non_test_files_in_test())

    if errors:
        for e in errors:
            print(e, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

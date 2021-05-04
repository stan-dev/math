#!/usr/bin/python

"""
Replacement for test-math-dependencies target in Makefile.

"""

from __future__ import print_function
import os
import sys
import re
import glob
from collections import defaultdict

winsfx = ".exe"
testsfx = "_test.cpp"
testsfx_no_ext = "_test"


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

def check_non_unique_test_names():
    test_files = files_in_folder("test/unit/")
    tests = {}
    duplicates = defaultdict(list)
    for filepath in test_files:
        if os.path.isfile(filepath) and filepath.endswith(testsfx):
            with open(filepath) as file:
                test_file_content = file.read()
                # look for TEST() and TEST_F()
                matches = re.findall(r"TEST(?:_F)?\((.*?)\)", test_file_content, re.DOTALL)
                for x in matches:
                    if not ("#" in x):
                        # strips for test names written in two lines
                        x_stripped = x.replace("\n", "").replace(" ", "").replace(",",", ")
                        if x_stripped in tests:                        
                            duplicates[x_stripped].append(filepath)
                        else:
                            tests[x_stripped] = filepath
    errors = []
    if len(duplicates)>0:
        duplicates_error_msg = ""
        for x in duplicates:
            duplicates_error_msg += "(" + x + ") in files:\n"
            # add the first found file first
            duplicates_error_msg += "\t" + tests[x] + "\n"
            for y in duplicates[x]:
                duplicates_error_msg += "\t" + y + "\n"
        errors.append(
            "Tests or test fixtures with non-unique names found in test/unit:\n\n" +
            duplicates_error_msg
        )
    return errors

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
    errors = []
    errors.extend(
        x + ":\n\t A .cpp file without the " + testsfx + " suffix found in test/unit/math/"
        for x in all_cpp
        if os.path.splitext(x)[1] == ".cpp" and x[-len(testsfx) :] != testsfx
    )
    errors.extend(
        x + ":\n\t A *_test.hpp file found. Must be *_test.cpp to run."
        for x in all_cpp
        if os.path.splitext(x)[1] == ".hpp" and x[-len(testsfx): -len(testsfx) + len(testsfx_no_ext)] == testsfx_no_ext
    )
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

    errors.extend(check_non_unique_test_names())

    if errors:
        for e in errors:
            print(e, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

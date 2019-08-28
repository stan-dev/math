#!/usr/bin/python

"""
Replacement for test-math-dependencies target in Makefile.

Call script with '-h' as an option to see a helpful message.
"""

from __future__ import print_function
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import os.path
import platform
import subprocess
import sys
import time
import glob
import re

winsfx = ".exe"
testsfx = "_test.cpp"

def grep_patterns(type, folder, patterns_and_messages, exclude_filters = []):
    """Checks the files in the provided folder for matches
    with any of the patterns. It returns an array of
    messages and the provided type with
    the line number. This check ignores comments.
    """
    errors = []
    folder.replace("/", os.sep)
    files = glob.glob(folder+os.sep+"**", recursive=True)
    filtered = []
    # filter out files that match exclude filters
    for exclude in exclude_filters:
        exclude.replace("/", os.sep)
        filtered += glob.glob(exclude, recursive=True)
    files = [x for x in files if x not in filtered]
    for f in files:
        if os.path.isfile(f):
            line_num = 0
            for line in open(f, "r"):
                line_num += 1
                # exclude comments
                # exclude line starting with "/*" or " * " and
                # if the matched patterns are behind "//"
                for p in patterns_and_messages:
                    if not re.search("^ \* |^/\*", line) and not re.search(".*//.*"+p["pattern"], line) and re.search(p["pattern"], line):
                            errors.append(f + " at line " + str(line_num) + ":\n\t" + "[" + type + "] " + p["message"])
    return errors

def check_non_test_files_in_test():
    all_cpp = glob.glob("test/unit/*/**", recursive=True)
    errors = ["Error: A .cpp file without the suffix " + testsfx + " in test/unit/*/*:\n\t" + x for x in all_cpp if x[-4:] == ".cpp" and x[-len(testsfx):] != testsfx]
    return errors

def main():
    errors = []
    # Check for files inside stan/math/prim that contain stan/math/rev or stan/math/fwd
    prim_checks = [
        { "pattern" : '<stan/math/rev/', "message" : 'File includes a stan/math/rev header file.'},
        { "pattern" : 'stan::math::var', "message" : 'File uses stan::math::var.'},
        { "pattern" : '<stan/math/fwd/', "message" : 'File includes a stan/math/fwd header file.'},
        { "pattern" : 'stan::math::fvar', "message" : 'File uses stan::math::fvar.'},
        { "pattern" : '<stan/math/mix/', "message" : 'File includes a stan/math/mix header file.'}
    ]
    errors += grep_patterns('prim', 'stan/math/prim', prim_checks)
    
    # Check for files inside stan/math/rev that contain stan/math/fwd or stan/math/mix
    rev_checks = [
        { "pattern" : '<stan/math/fwd/', "message" : 'File includes a stan/math/fwd header file.'},
        { "pattern" : 'stan::math::fvar', "message" : 'File uses stan::math::fvar.'},
        { "pattern" : '<stan/math/mix/', "message" : 'File includes a stan/math/mix header file.'}
    ]
    errors += grep_patterns('rev', 'stan/math/rev', rev_checks)

    # Check for files inside stan/math/*/scal that contain stan/math/*/arr or stan/math/*/mat
    scal_checks = [
        { "pattern" : '<stan/math/.*/arr/', "message" : 'File includes an array header file.'},
        { "pattern" : '<vector>', "message" : 'File includes an std::vector header.'},
        { "pattern" : 'std::vector', "message" : 'File uses std::vector.'}, 
        { "pattern" : '<stan/math/.*/mat/', "message" : 'File includes a matrix header file.'}, 
        { "pattern" : '<Eigen', "message" : 'File includes an Eigen header.'},
        { "pattern" : 'Eigen::', "message" : 'File uses Eigen.'}
    ]
    errors += grep_patterns('scal', 'stan/math/*/scal', scal_checks)

    # Check for files inside stan/math/*/arr that contain stan/math/*/mat or Eigen
    arr_checks = [
        { "pattern" : '<stan/math/.*/mat/', "message" : 'File includes an matrix header file.'},
        { "pattern" : '<Eigen', "message" : 'File includes an Eigen header.'},
        { "pattern" : 'Eigen::', "message" : 'File uses Eigen.'}
    ]    
    errors += grep_patterns('arr', 'stan/math/*/arr', arr_checks)

    # Check to make sure we use C++14 constructs in stan/math
    cpp14_checks = [
        { "pattern" : 'boost::is_unsigned<', "message" : 'File uses boost::is_unsigned instead of std::is_unsigned.'},
        { "pattern" : '<boost/type_traits/is_unsigned>',\
            "message" : 'File includes <boost/type_traits/is_unsigned.hpp> instead of <type_traits>.'},
        { "pattern" : 'boost::is_arithmetic<', "message" : 'File uses boost::is_arithmetic instead of std::is_arithmetic.'},
        { "pattern" : '<boost/type_traits/is_arithmetic.hpp>',\
            "message" : 'File includes <boost/type_traits/is_unsigned.hpp> instead of <type_traits>.'},
        { "pattern" : 'boost::is_convertible<', "message" : 'File uses boost::is_convertible instead of std::is_convertible.'},
        { "pattern" : '<boost/type_traits/is_convertible.hpp>',\
            "message" : 'File includes <boost/type_traits/is_convertible.hpp> instead of <type_traits>.'},
        { "pattern" : 'boost::is_same<', "message" : 'File uses boost::is_same instead of std::is_same.'},
        { "pattern" : '<boost/type_traits/is_same.hpp>',\
            "message" : 'File includes <boost/type_traits/is_same.hpp> instead of <type_traits>.'},
        { "pattern" : 'boost::enable_if_c<', "message" : 'File uses boost::enable_if_c instead of std::enable_if.'},
        { "pattern" : 'boost::enable_if<', "message" : 'File uses boost::enable_if instead of std::enable_if.'},
        { "pattern" : 'boost::disable_if<', "message" : 'File uses boost::disable_if instead of std::enable_if.'},
        { "pattern" : '<boost/utility/enable_if.hpp>',\
            "message" : 'Replace \<boost/utility/enable_if.hpp\> with \<type_traits\>.'}
    ]    
    errors += grep_patterns('C++14', 'stan/math', cpp14_checks)
    
    # Check for includes of stan/math/*/meta/*.hpp inside stan/math, excluding meta.hpp files and the /meta subfolder
    meta_checks = [
        { "pattern" : '<stan/math/.*/meta/.*hpp', "message" : 'File includes */meta/*.hpp header file. Should include meta.hpp'}
    ]
    meta_exclude = ['stan/math/*/*/meta/*', 'stan/math/*/meta.hpp']
    errors += grep_patterns('meta', 'stan/math', meta_checks, meta_exclude)

    errors += check_non_test_files_in_test()

    if(len(errors) > 0):
        for e in errors:
            print(e)
        sys.exit(1)
        
if __name__ == "__main__":
    main()

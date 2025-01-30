"""

    Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.
    Copyright (c) Sorin Fetche

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

    This program generates a single header file from a file including multiple C/C++ headers.

    Usage:

        python3 generate_single_header.py  --help

        e.g. python3 generate_single_header.py -i include/boost/leaf/detail/all.hpp -p include  -o test/leaf.hpp boost/leaf

"""

import argparse
import os
import filecmp
import re
from datetime import date
import subprocess

def compare_and_update(old_file, new_file):
    if not os.path.exists(old_file):
        os.rename(new_file, old_file)
    elif filecmp.cmp(old_file, new_file, shallow=False):
        os.remove(new_file)
    else:
        os.remove(old_file)
        os.rename(new_file, old_file)

included = {}
total_line_count = 14
def append(input_file_name, input_file, output_file, regex_includes, include_folder, line_directive_prefix, depth):
    global total_line_count
    inside_copyright = False
    line_count = 1
    for line in input_file:
        line_count += 1
        if 'Emil Dotchevski' in line:
            inside_copyright = True
        if inside_copyright:
            if line.startswith('//') :
                continue
            if not line.strip():
                inside_copyright = False
                if depth > 0:
                    output_file.write('%s#line %d "%s"\n' % (line_directive_prefix, line_count, input_file_name))
                continue
        result = regex_includes.search(line)
        if result:
            next_input_file_name = result.group('include')
            if next_input_file_name in included:
                output_file.write('// %s // Expanded at line %d\n' % (line.strip(), included[next_input_file_name]))
                total_line_count += 1
            else:
                included[next_input_file_name] = total_line_count
                print('%s (%d)' % (next_input_file_name, total_line_count))
                with open(os.path.join(include_folder, next_input_file_name), 'r') as next_input_file:
                    output_file.write('// >>> %s' % (line))
                    total_line_count += 2
                    append(next_input_file_name, next_input_file, output_file, regex_includes, include_folder, line_directive_prefix, depth + 1)
                    if depth > 0:
                        output_file.write('// <<< %s%s#line %d "%s"\n' % (line, line_directive_prefix, line_count, input_file_name))
                        total_line_count += 2
        else:
            output_file.write(line)
            total_line_count += 1

def _main():
    parser = argparse.ArgumentParser(
        description='Generates a single include file from a file including multiple C/C++ headers')
    parser.add_argument('-i', '--input', action='store', type=str, default='in.cpp',
        help='Input file including the headers to be merged')
    parser.add_argument('-o', '--output', action='store', type=str, default='out.cpp',
        help='Output file. NOTE: It will be overwritten!')
    parser.add_argument('-p', '--path', action='store', type=str, default='.',
        help='Include path')
    parser.add_argument('--hash', action='store', type=str,
        help='The git hash to print in the output file, e.g. the output of "git rev-parse HEAD"')
    parser.add_argument('prefix', action='store', type=str,
        help='Non-empty include file prefix (e.g. a/b)')
    parser.add_argument('--linerefs', action='store_true',
        help='Output #line references to the original files. By default the line references are written commented out')
    args = parser.parse_args()

    regex_includes = re.compile(r"""^\s*#[\t\s]*include[\t\s]*("|\<)(?P<include>%s.*)("|\>)""" % args.prefix)
    print('Rebuilding %s:' % args.input)
    tmp_file_name = args.output + '.tmp'
    with open(tmp_file_name, 'w') as tmp_file, open(args.input, 'r') as input_file:
        tmp_file.write(
            '#ifndef BOOST_LEAF_HPP_INCLUDED\n'
            '#define BOOST_LEAF_HPP_INCLUDED\n'
            '\n'
            '// Boost LEAF single header distribution. Do not edit.\n'
            '// Generated on ' + date.today().strftime('%b %d, %Y'))
        if args.hash:
            tmp_file.write(
                ' from https://github.com/boostorg/leaf/tree/' + args.hash[0:7])
        tmp_file.write(
            '.\n'
            '\n'
            '// Latest published version of this file: https://raw.githubusercontent.com/boostorg/leaf/gh-pages/leaf.hpp.\n'
            '\n'
            '// Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.\n'
            '// Distributed under the Boost Software License, Version 1.0. (See accompanying\n'
            '// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)\n'
            '\n')
        append(args.input, input_file, tmp_file, regex_includes, args.path, '' if args.linerefs else '// ', 0)
        tmp_file.write(
            '\n'
            '#endif // BOOST_LEAF_HPP_INCLUDED\n' )
    compare_and_update(args.output, tmp_file_name)
#       print('%d' % total_line_count)

if __name__ == "__main__":
    _main()

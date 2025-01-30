#!/usr/bin/env python3

# Copyright (c) 2024 T. Zachary Laine
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Get the latest version of the data file at:
# https://www.unicode.org/Public/UCD/latest/ucd/CaseFolding.txt

import itertools

f = open('CaseFolding.txt')

lines = f.readlines()

lines = filter(
    lambda x: False if (x.startswith('#') or x == '\n' or ' T;' in x or ' S;' in x) else True,
    lines)

lines = list(map(lambda x: x.split(';')[:-1], lines))

f_line_indices = list(filter(lambda i: ' F;' in lines[i], range(len(lines))))

# Remove all C; lines that come before a F; line that applies to the same code
# point.
for f_idx in reversed(f_line_indices):
    prev = f_idx - 1
    if 0 <= prev:
        if lines[f_idx][0] == lines[prev][0]:
            lines = lines[:prev] + lines[prev + 1:]

max_mapping_len = 1
for line in lines:
    mappings = line[2].strip().split(' ')
    max_mapping_len = max(len(mappings), max_mapping_len)

def print_cps():
    print('// Initial code points from CaseFolding.txt')
    print('char32_t const cps[] = {')
    str_ = ''
    idx = 0
    max_per_line = 8
    for cp in map(lambda x: '0x' + x[0], lines):
        str_ += cp + ', '
        idx += 1
        if (idx % max_per_line) == 0:
            print('    ' + str_)
            str_ = ''
    if idx % max_per_line:
        print('    ' + str_)
    print('};\n')
    print(f'[[maybe_unused]] char32_t const max_test_cp = 0x{lines[-1][0]} + 100;\n')


array_t = f'std::array<uint32_t, {max_mapping_len} + 1>'

def print_test(line):
    num_mapping_cps = len(line[2].strip().split(' '))
    mapping = ', '.join(map(lambda x: '0x' + x, line[2].strip().split(' ')))
    trailing = ', 0' * (max_mapping_len + 1 - num_mapping_cps)
    zeros = ', 0' * (max_mapping_len)
    print(f'''    {{
        {array_t} const expected = {{ {mapping}{trailing} }};
        {array_t} result = {{ 0{zeros} }};
        boost::parser::detail::case_fold(0x{line[0]}, result.begin());
        BOOST_TEST(result == expected);
    }}''')

def print_tests():
    idx = 0
    max_per_TEST = 50
    print(f'// hits_{int(idx / max_per_TEST)})\n{{')
    for line in lines:
        idx += 1
        if (idx % max_per_TEST) == 0:
            print(f'''}}

// test_{int(idx / max_per_TEST)})
{{''')
        print_test(line)
    print('}\n')

print('''// Copyright (c) 2024 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Warning: This header is auto-generated (see misc/generate_case_fold_tests.py).

#include <boost/parser/parser.hpp>

#include <boost/core/lightweight_test.hpp>


int main()
{
''')

print_cps()
print_tests()

print(f'''// misses
{{
    char32_t next_cp = 0;
    for (char32_t const * it = cps; it != std::end(cps); ++it) {{
        for (char32_t cp = next_cp; cp < *it; ++cp) {{
            {array_t} const expected = {{ cp, 0 }};
            {array_t} result = {{ 0 }};
            auto const first = result.data();
            auto const last = boost::parser::detail::case_fold(cp, first);
            BOOST_TEST(std::equal(first, last, expected.begin()));
            BOOST_TEST(result == expected);
        }}
        next_cp = *it + 1;
    }}
}}

return boost::report_errors();
}}''')

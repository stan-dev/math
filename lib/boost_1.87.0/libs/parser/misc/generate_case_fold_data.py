#!/usr/bin/env python3

# Copyright (c) 2024 T. Zachary Laine
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Get the latest version of the data file at:
# https://www.unicode.org/Public/UCD/latest/ucd/CaseFolding.txt

import itertools

def print_structs():
    print('''struct short_mapping_range {
    char32_t cp_first_;
    char32_t cp_last_;
    uint16_t stride_;
    uint16_t first_idx_;
};''')

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

# The mapping code points, in the order listed in lines (only individual mapping-cps).
all_mapping_cps = []
max_mapping_len = 1
# The mapping code points, in the order listed in lines (only multiple mapping-cps).
big_mappings = []
for line in lines:
    mappings = line[-1].strip().split(' ')
    max_mapping_len = max(len(mappings), max_mapping_len)
    if 1 < len(mappings):
        big_mappings.append((line[0], mappings))
    else:
        line.append(len(all_mapping_cps))
        all_mapping_cps += mappings

def print_long_mapping():
    print(f'''inline constexpr int longest_mapping = {max_mapping_len};
struct long_mapping {{
    char32_t cp_;
    char32_t mapping_[longest_mapping + 1];
}};''')
    print(f'inline constexpr long_mapping  long_mappings[{len(big_mappings)}] = {{')
    for big_mapping in big_mappings:
        mappings = ', '.join(map(lambda x: '0x' + x, big_mapping[1]))
        print(f'    {{0x{big_mapping[0]}, {{{mappings} , 0}}}},')
    print('};')

def print_mapping_cps():
    print(f'inline constexpr char32_t single_mapping_cps[{len(all_mapping_cps)}] = {{')
    str_ = ''
    idx = 0
    max_per_line = 8
    for cp in map(lambda x: '0x' + x, all_mapping_cps):
        str_ += cp + ', '
        idx += 1
        if (idx % max_per_line) == 0:
            print('    ' + str_)
            str_ = ''
    if idx % max_per_line:
        print('    ' + str_)
    print('};\n')

# Each range consists of a first-index, last-index, and the stride of values
# within [first, last).
ranges = []

lines = list(filter(lambda x: len(x[2].strip().split(' ')) == 1, lines))
pairs = zip(lines[:-1], lines[1:])

def pair_diff(pair):
    prev_cp = int('0x' + pair[0][0], 0)
    cp = int('0x' + pair[1][0], 0)
    return cp - prev_cp

for k, g in itertools.groupby(pairs, pair_diff):
    g = list(g)
    ranges.append((g[0][0][0], g[-1][1][0], k, g[0][0][-1]))
ranges[-1] = (ranges[-1][0], lines[-1][0] + ' + 1', ranges[-1][2], ranges[-1][3])

def print_ranges():
    print(f'inline constexpr short_mapping_range mapping_ranges[{len(ranges)}] = {{')
    for r in ranges:
        print(f'    {{0x{r[0]}, 0x{r[1]}, {r[2]}, {r[3]}}}, ')
    print('};')

print('''// Copyright (c) 2024 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Warning: This header is auto-generated (see misc/generate_case_fold_data.py).
// The lack of include guards is intentional.

namespace boost::parser::detail {

''')

print_structs()
print_long_mapping()
print_ranges()
print_mapping_cps()

print('}')

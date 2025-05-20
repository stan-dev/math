#
# Copyright (c) 2023 Alan de Freitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/url
#

import re


def replace_snippets(input_file, output_file):
    with open(input_file, 'r') as input_f:
        code = input_f.read()

    # Define a regular expression pattern to match snippet headers
    snippet_pattern = r'\/\/\[\s*(\w+)\s*'
    current_snippet = None
    snippet_close_pattern = r'\/\/\]'
    output_lines = []
    for line in code.splitlines():
        match = re.search(snippet_pattern, line)
        if match:
            output_lines.append(re.sub(snippet_pattern, r'// tag::\1[]', line))
            current_snippet = match.group(1)
            continue
        # Check if the line matches the closing snippet pattern
        match = re.search(snippet_close_pattern, line)
        if match:
            output_lines.append(re.sub(snippet_close_pattern, f'// end::{current_snippet}[]', line))
            continue
        output_lines.append(line)

    with open(output_file, 'w') as output_f:
        output_f.write('\n'.join(output_lines))


if __name__ == "__main__":
    input_file = 'snippets.cpp'
    output_file = 'snippets_out.cpp'
    replace_snippets(input_file, output_file)

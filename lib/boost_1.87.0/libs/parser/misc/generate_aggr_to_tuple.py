#!/usr/bin/python3

# Copyright (c) 2023 T. Zachary Laine
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

print('''// Copyright (c) 2023 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Warning: This header is auto-generated (see misc/generate_aggr_to_tuple.py).
// The lack of include guards is intentional.

namespace boost::parser::detail {

''')

specializations = 100
bindings = []

def bindings_str(fn, max_bindings_on_a_line):
    num_bindings = len(bindings)
    mapped = map(fn, bindings)
    retval = ''
    for i in range(int(num_bindings / max_bindings_on_a_line) + 1):
        lo = i * max_bindings_on_a_line
        if num_bindings <= lo:
            return retval
        hi = min(lo + max_bindings_on_a_line, num_bindings)
        retval += '    ' + ', '.join(map(fn, bindings[lo:hi]))
        if hi < num_bindings:
            retval += ','
        retval += '\n'
    return retval

def destructuring_str():
    return f'''auto & [
{bindings_str(lambda x: x, 15)}] = x;'''

def typenames_str():
    return f'''<
{bindings_str(lambda x: f'decltype({x}) &', 5)}>'''

def initializers_str():
    return f'''return parser::tuple{typenames_str()}(
{bindings_str(lambda x: f'{x}', 15)});'''

for i in range(1, specializations + 1):
    bindings.append('_{:02x}'.format(i))

    print(f'''template<> struct tie_aggregate_impl<{i}> {{
template<typename T> static constexpr auto call(T & x) {{''')
    print(destructuring_str())
    print(initializers_str())
    print('''}
};

''')

print('}\n')

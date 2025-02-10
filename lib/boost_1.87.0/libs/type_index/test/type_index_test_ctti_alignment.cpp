//
// Copyright 2012-2024 Antony Polukhin.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/type_index/ctti_type_index.hpp>

int main() {
    static_assert(
        alignof(boost::typeindex::detail::ctti_data) == alignof(char),
        "Alignments of boost::typeindex::detail::ctti_data and char differ. "
        "It is unsafe to reinterpret_cast between them."
    );
}


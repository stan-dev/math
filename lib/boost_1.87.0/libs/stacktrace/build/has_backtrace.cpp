// Copyright Antony Polukhin, 2016-2024.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifdef BOOST_STACKTRACE_BACKTRACE_INCLUDE_FILE
#   include BOOST_STACKTRACE_BACKTRACE_INCLUDE_FILE
#else
#   include <backtrace.h>
#endif

#include <unwind.h>

int main() {
    backtrace_state* state = backtrace_create_state(
        0, 1, 0, 0
    );
    (void)state;
}

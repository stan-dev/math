//  Copyright Antony Polukhin, 2023-2024.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).

#include <boost/lexical_cast.hpp>

#include <boost/core/lightweight_test.hpp>

#include <boost/lexical_cast.hpp>

struct oops {
    operator int () const {
        return 41;
    }
};

inline std::ostream& operator<<(std::ostream& os, const oops&) {
    return os << 42;
}

int main () {
    auto val = boost::lexical_cast<int>(oops{});
    BOOST_TEST_EQ(val, 42);

    return boost::report_errors();
}

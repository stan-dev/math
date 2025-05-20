//  Unit test for boost::lexical_cast.
//
//  See http://www.boost.org for most recent version, including documentation.
//
//  Copyright Alexander Nasonov, 2007.
//  Copyright Antony Polukhin, 2023-2024.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).
//
//  Test that Source can be non-copyable.

#include <boost/lexical_cast.hpp>
#include <boost/core/noncopyable.hpp>

#include <boost/core/lightweight_test.hpp>

class Noncopyable : private boost::noncopyable
{
public:
    Noncopyable() {}
};

inline std::ostream &operator<<(std::ostream &out, const Noncopyable&)
{
    return out << "Noncopyable";
}

void test_noncopyable()
{
    Noncopyable x;
    BOOST_TEST(boost::lexical_cast<std::string>(x) == "Noncopyable");
}

int main()
{
    test_noncopyable();
    return boost::report_errors();
}

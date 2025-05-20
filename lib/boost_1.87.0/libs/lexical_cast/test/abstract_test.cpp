//  Unit test for boost::lexical_cast.
//
//  See http://www.boost.org for most recent version, including documentation.
//
//  Copyright Sergey Shandar 2005, Alexander Nasonov, 2007.
//  Copyright Antony Polukhin, 2023-2024.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).
//
// Test abstract class. Bug 1358600:
// http://sf.net/tracker/?func=detail&aid=1358600&group_id=7586&atid=107586

#include <boost/lexical_cast.hpp>

#include <boost/core/lightweight_test.hpp>

class A
{
public:
    virtual void out(std::ostream &) const = 0;
    virtual ~A() {}
};

class B: public A
{
public:
    virtual void out(std::ostream &O) const { O << "B"; }
};

std::ostream &operator<<(std::ostream &O, const A &a)
{
    a.out(O);
    return O;
}

void test_abstract()
{
    const A &a = B();
    BOOST_TEST(boost::lexical_cast<std::string>(a) == "B");
}

int main()
{
    test_abstract();
    return boost::report_errors();
}


//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

// Test that header file is self-contained.
#include <boost/json/visit.hpp>

#include "test_suite.hpp"

namespace boost {
namespace json {

class visit_test
{
public:
    template<class T>
    void
    check_const(kind k, T t)
    {
        struct f
        {
            kind k;
            bool operator()(std::nullptr_t const&) const { return k == kind::null; }
            bool operator()(bool const&) const { return k == kind::bool_; }
            bool operator()(std::int64_t const&) const { return k == kind::int64; }
            bool operator()(std::uint64_t const&) const { return k == kind::uint64; }
            bool operator()(double const&) const { return k == kind::double_; }
            bool operator()(string const&) const { return k == kind::string; }
            bool operator()(array const&) const { return k == kind::array; }
            bool operator()(object const&) const { return k == kind::object; }
            bool operator()(...) const { return false; }
        };
        value const v(t);
        BOOST_TEST(visit(f{k}, v));
    }

    template<class T>
    void
    check_mutable(kind k, T t)
    {
        struct f
        {
            kind k;
            bool operator()(std::nullptr_t&) { return k == kind::null; }
            bool operator()(bool&) { return k == kind::bool_; }
            bool operator()(std::int64_t&) { return k == kind::int64; }
            bool operator()(std::uint64_t&) { return k == kind::uint64; }
            bool operator()(double&) { return k == kind::double_; }
            bool operator()(string&) { return k == kind::string; }
            bool operator()(array&) { return k == kind::array; }
            bool operator()(object&) { return k == kind::object; }
            bool operator()(...) { return false; }
        };
        value v(t);
        BOOST_TEST(visit(f{k}, v));
    }

    template<class T>
    void
    check_rvalue(kind k, T t)
    {
        struct f
        {
            kind k;
            bool operator()(std::nullptr_t&&) && { return k == kind::null; }
            bool operator()(bool&&) && { return k == kind::bool_; }
            bool operator()(std::int64_t&&) && { return k == kind::int64; }
            bool operator()(std::uint64_t&&) && { return k == kind::uint64; }
            bool operator()(double&&) && { return k == kind::double_; }
            bool operator()(string&&) && { return k == kind::string; }
            bool operator()(array&&) && { return k == kind::array; }
            bool operator()(object&&) && { return k == kind::object; }
            bool operator()(...) && { return false; }
        };
        value v(t);
        f f_{k};
        BOOST_TEST(visit( std::move(f_), std::move(v) ));
    }

    template<class T>
    void
    check_nonref(kind k, T t)
    {
        struct f
        {
            kind k;
            bool operator()(std::nullptr_t) && { return k == kind::null; }
            bool operator()(bool) && { return k == kind::bool_; }
            bool operator()(std::int64_t) && { return k == kind::int64; }
            bool operator()(std::uint64_t) && { return k == kind::uint64; }
            bool operator()(double) && { return k == kind::double_; }
            bool operator()(string) && { return k == kind::string; }
            bool operator()(array) && { return k == kind::array; }
            bool operator()(object) && { return k == kind::object; }
            bool operator()(...) && { return false; }
        };
        value v(t);
        f f_{k};
        BOOST_TEST(visit( std::move(f_), std::move(v) ));
    }

    void
    testVisit()
    {
        check_const(kind::null,    nullptr);
        check_const(kind::bool_,   true);
        check_const(kind::int64,   -1);
        check_const(kind::uint64,  1U);
        check_const(kind::double_, 3.14);
        check_const(kind::string,  string_kind);
        check_const(kind::array,   array_kind);
        check_const(kind::object,  object_kind);

        check_mutable(kind::null,    nullptr);
        check_mutable(kind::bool_,   true);
        check_mutable(kind::int64,   -1);
        check_mutable(kind::uint64,  1U);
        check_mutable(kind::double_, 3.14);
        check_mutable(kind::string,  string_kind);
        check_mutable(kind::array,   array_kind);
        check_mutable(kind::object,  object_kind);

        check_rvalue(kind::null,    nullptr);
        check_rvalue(kind::bool_,   true);
        check_rvalue(kind::int64,   -1);
        check_rvalue(kind::uint64,  1U);
        check_rvalue(kind::double_, 3.14);
        check_rvalue(kind::string,  string_kind);
        check_rvalue(kind::array,   array_kind);
        check_rvalue(kind::object,  object_kind);

        check_nonref(kind::null,    nullptr);
        check_nonref(kind::bool_,   true);
        check_nonref(kind::int64,   -1);
        check_nonref(kind::uint64,  1U);
        check_nonref(kind::double_, 3.14);
        check_nonref(kind::string,  string_kind);
        check_nonref(kind::array,   array_kind);
        check_nonref(kind::object,  object_kind);
    }

    void run()
    {
        testVisit();
    }
};

TEST_SUITE(visit_test, "boost.json.visit");

} // namespace json
} // namespace boost

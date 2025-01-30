// Copyright (c) 2024 Mikhail Khachayants (mkhachaiants@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

#include <boost/json.hpp>

#if !defined(BOOST_DESCRIBE_CXX14)

#include <boost/config/pragma_message.hpp>

BOOST_PRAGMA_MESSAGE( "This example requires C++14" )

int main() {}

#else

#include <boost/json/parse_into.hpp>
#include <boost/variant2/variant.hpp>
#include <boost/describe.hpp>
#include <map>

#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
# include <optional>
# define IF_CXX17_HDR_OPTIONAL(...) __VA_ARGS__
#else
# define IF_CXX17_HDR_OPTIONAL(...)
#endif // BOOST_NO_CXX17_HDR_OPTIONAL

using namespace boost::json;

struct Object
{
    bool b;
    float f;
    double d;
    std::int64_t i64;
    std::uint64_t u64;
    std::string s;
    std::vector<bool> v1;
    std::vector<std::int64_t> v2;
    std::vector<std::uint64_t> v3;
    std::array<bool, 3> a1;
    std::array<std::int64_t, 3> a2;
    std::array<std::uint64_t, 3> a3;
    std::map<std::string, std::int64_t> m1;
    std::map<std::string, std::string> m2;
    std::map<std::string, double> m3;
    std::tuple<bool, std::uint64_t, std::int64_t, double, std::string> t1;
    std::tuple<std::array<std::string, 3>, std::array<double, 3>, std::nullptr_t> t2;
    std::tuple<std::vector<std::string>, std::vector<double>> t3;
    boost::variant2::variant<bool, std::uint64_t, std::int64_t, double, std::string> v;

#ifndef BOOST_NO_CXX17_HDR_OPTIONAL
    std::optional<bool> ob;
    std::optional<std::int64_t> oi;
    std::optional<std::uint64_t> ou;
    std::optional<double> od;
    std::optional<std::string> os;
#endif // BOOST_NO_CXX17_HDR_OPTIONAL
};

BOOST_DESCRIBE_STRUCT(Object, (),
    (b, i64, u64, f, d, s, v1, v2, v3, a1, a2, a3, m1, m2, m3, t1, t2, t3, v,
    IF_CXX17_HDR_OPTIONAL(ob, oi, ou, od, os)))


bool
fuzz_direct_parse(string_view sv)
{
    Object object;
    boost::system::error_code ec;
    parse_into(object, sv, ec);
    return !ec;
}

extern "C"
int
LLVMFuzzerTestOneInput(
        const uint8_t* data, size_t size)
{
    try
    {
        string_view sv{reinterpret_cast<
            const char*>(data), size};
        fuzz_direct_parse(sv);
    }
    catch(...)
    {
    }
    return 0;
}

#endif // !defined(BOOST_DESCRIBE_CXX14)

/**
 *   Copyright (C) 2018 T. Zachary Laine
 *
 *   Distributed under the Boost Software License, Version 1.0. (See
 *   accompanying file LICENSE_1_0.txt or copy at
 *   http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/parser/tuple.hpp>

#include <boost/core/lightweight_test.hpp>

#include "ill_formed.hpp"

#include <string>


namespace bp = boost::parser;


struct empty
{};

struct int_
{
    int i;
};

struct int_string
{
    int i;
    std::string s;
};

struct string_int
{
    std::string s;
    int i;
};

struct string_int_double
{
    std::string s;
    int i;
    double d;
};

struct string_int_double_priv
{
    std::string s;

private:
    [[maybe_unused]] int i;
    [[maybe_unused]] double d;
};

struct string_int_double_no_copy_move
{
    string_int_double_no_copy_move(string_int_double_no_copy_move const &) =
        delete;
    string_int_double_no_copy_move(string_int_double_no_copy_move &&) = delete;

private:
    std::string s;
    [[maybe_unused]] int i;
    [[maybe_unused]] double d;
};

struct employee
{
    int age;
    std::string surname;
    std::string forename;
    double salary;
};

int main()
{

// struct_arity
{
    using bp::detail::struct_arity_v;

    static_assert(struct_arity_v<empty> == 0);
    static_assert(struct_arity_v<int_> == 1);
    static_assert(struct_arity_v<int_string> == 2);
    static_assert(struct_arity_v<string_int> == 2);
    static_assert(struct_arity_v<string_int_double> == 3);

    // 1 due to the copy ctor.
    static_assert(struct_arity_v<string_int_double_priv> == 1);

    static_assert(struct_arity_v<string_int_double_no_copy_move> <= 0);
}

// is_struct_assignable
{
    using bp::detail::is_struct_assignable_v;

    static_assert(!is_struct_assignable_v<void, void>);
    static_assert(!is_struct_assignable_v<int, bp::tuple<int, std::string>>);
    static_assert(!is_struct_assignable_v<void, bp::tuple<int, std::string>>);
    static_assert(!is_struct_assignable_v<empty, bp::tuple<int, std::string>>);
    static_assert(!is_struct_assignable_v<int_, bp::tuple<int, std::string>>);
    static_assert(
        !is_struct_assignable_v<string_int, bp::tuple<int, std::string>>);
    static_assert(!is_struct_assignable_v<
                  string_int_double,
                  bp::tuple<std::string, int>>);
    static_assert(!is_struct_assignable_v<
                  string_int,
                  bp::tuple<std::string, int, double>>);
    static_assert(!is_struct_assignable_v<
                  string_int_double_priv,
                  bp::tuple<std::string, int, double>>);
    static_assert(!is_struct_assignable_v<
                  string_int_double_no_copy_move,
                  bp::tuple<std::string, int, double>>);

    static_assert(
        is_struct_assignable_v<int_string, bp::tuple<int, std::string>>);
    static_assert(
        is_struct_assignable_v<int_string, bp::tuple<short, char const *>>);

    static_assert(
        is_struct_assignable_v<string_int, bp::tuple<std::string, int>>);
    static_assert(
        is_struct_assignable_v<string_int, bp::tuple<char const *, short>>);

    static_assert(is_struct_assignable_v<
                  string_int_double,
                  bp::tuple<std::string, int, double>>);
    static_assert(is_struct_assignable_v<
                  string_int_double,
                  bp::tuple<char const *, short, float>>);
}

// is_tuple_assignable
{
    using bp::detail::is_tuple_assignable_v;

    static_assert(!is_tuple_assignable_v<void, void>);
    static_assert(!is_tuple_assignable_v<bp::tuple<int, std::string>, int>);
    static_assert(!is_tuple_assignable_v<bp::tuple<int, std::string>, void>);
    static_assert(!is_tuple_assignable_v<bp::tuple<int, std::string>, empty>);
    static_assert(!is_tuple_assignable_v<bp::tuple<int, std::string>, int_>);
    static_assert(
        !is_tuple_assignable_v<bp::tuple<int, std::string>, string_int>);
    static_assert(
        !is_tuple_assignable_v<bp::tuple<std::string, int>, string_int_double>);
    static_assert(!is_tuple_assignable_v<
                  bp::tuple<std::string, int, double>,
                  string_int>);
    static_assert(!is_tuple_assignable_v<
                  bp::tuple<std::string, int, double>,
                  string_int_double_priv>);
    static_assert(!is_tuple_assignable_v<
                  bp::tuple<std::string, int, double>,
                  string_int_double_no_copy_move>);

    static_assert(
        std::is_aggregate_v<int_string> &&
        bp::detail::struct_arity_v<int_string> ==
            bp::detail::tuple_size_<bp::tuple<int, std::string>>);

    static_assert(
        is_tuple_assignable_v<bp::tuple<int, std::string>, int_string>);
    static_assert(
        is_tuple_assignable_v<bp::tuple<short, std::string>, int_string>);

    static_assert(
        is_tuple_assignable_v<bp::tuple<std::string, int>, string_int>);
    static_assert(
        is_tuple_assignable_v<bp::tuple<std::string, short>, string_int>);

    static_assert(is_tuple_assignable_v<
                  bp::tuple<std::string, int, double>,
                  string_int_double>);
    static_assert(is_tuple_assignable_v<
                  bp::tuple<std::string, short, float>,
                  string_int_double>);
}

// tuple_to_aggregate
{
    {
        employee const expected_employee = {32, "Last", "First", 50000.0};
        bp::tuple<int, std::string, std::string, double> tup(
            32, "Last", "First", 50000.0);

        employee assignee = bp::detail::tuple_to_aggregate<employee>(
            std::move(tup),
            std::make_integer_sequence<
                int,
                bp::detail::tuple_size_<decltype(tup)>>());

        BOOST_TEST(assignee.age == expected_employee.age);
        BOOST_TEST(assignee.surname == expected_employee.surname);
        BOOST_TEST(assignee.forename == expected_employee.forename);
        BOOST_TEST(assignee.salary == expected_employee.salary);

        // Verify move.
        BOOST_TEST(bp::get(tup, bp::llong<1>{}).empty());
        BOOST_TEST(bp::get(tup, bp::llong<2>{}).empty());
    }
    {
        employee const expected_employee = {32, "Last", "First", 50000.0};
        bp::tuple<int, char const *, char const *, double> tup(
            32, "Last", "First", 50000.0);

        employee assignee = bp::detail::tuple_to_aggregate<employee>(
            std::move(tup),
            std::make_integer_sequence<
                int,
                bp::detail::tuple_size_<decltype(tup)>>());

        BOOST_TEST(assignee.age == expected_employee.age);
        BOOST_TEST(assignee.surname == expected_employee.surname);
        BOOST_TEST(assignee.forename == expected_employee.forename);
        BOOST_TEST(assignee.salary == expected_employee.salary);
    }
}

// tie_aggregate
{
    employee assigner = {32, "Last", "First", 50000.0};

    auto tup = bp::detail::tie_aggregate(assigner);

    BOOST_TEST(bp::get(tup, bp::llong<0>{}) == 32);
    BOOST_TEST(bp::get(tup, bp::llong<1>{}) == "Last");
    BOOST_TEST(bp::get(tup, bp::llong<2>{}) == "First");
    BOOST_TEST(bp::get(tup, bp::llong<3>{}) == 50000.0);
}

return boost::report_errors();
}

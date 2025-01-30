//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/field_view.hpp>
#include <boost/mysql/rows.hpp>
#include <boost/mysql/rows_view.hpp>

#include <boost/mysql/detail/access.hpp>
#include <boost/mysql/detail/rows_iterator.hpp>

#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <tuple>
#include <utility>

#include "test_common/create_basic.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;

BOOST_TEST_DONT_PRINT_LOG_VALUE(boost::mysql::detail::rows_iterator)

BOOST_AUTO_TEST_SUITE(test_rows_iterator)

// Provide a uniform interface for both rows and rows_view types,
// so we can use template tests to reduce duplication
struct rows_view_wrapper
{
    std::vector<field_view> fields;
    rows_view r;

    rows_view_wrapper() = default;

    template <typename... Args>
    rows_view_wrapper(std::size_t num_columns, Args&&... args)
        : fields(make_fv_vector(std::forward<Args>(args)...)),
          r(makerowsv(fields.data(), sizeof...(args), num_columns))
    {
    }

    using iterator = rows_view::const_iterator;

    rows_view::const_iterator begin() const noexcept { return r.begin(); }
    rows_view::const_iterator end() const noexcept { return r.end(); }
};

struct rows_wrapper
{
    rows r;

    rows_wrapper() = default;

    template <typename... Args>
    rows_wrapper(std::size_t num_columns, Args&&... args)
    {
        auto fields = make_fv_arr(std::forward<Args>(args)...);
        r = makerowsv(fields.data(), fields.size(), num_columns);
    }

    using iterator = rows::const_iterator;

    rows::const_iterator begin() const noexcept { return r.begin(); }
    rows::const_iterator end() const noexcept { return r.end(); }
};

using rows_types = std::tuple<rows_view_wrapper, rows_wrapper>;

BOOST_AUTO_TEST_SUITE(range_iteration)
BOOST_AUTO_TEST_CASE_TEMPLATE(empty, RowType, rows_types)
{
    RowType r;
    std::vector<row_view> fv(r.begin(), r.end());
    BOOST_TEST(fv.size() == 0u);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(one_row_one_column, RowType, rows_types)
{
    RowType r(1, 42);
    std::vector<row_view> fv(r.begin(), r.end());
    BOOST_TEST(fv.size() == 1u);
    BOOST_TEST(fv[0] == makerow(42));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(one_row_several_columns, RowType, rows_types)
{
    RowType r(2, 80u, "abc");
    std::vector<row_view> fv(r.begin(), r.end());
    BOOST_TEST(fv.size() == 1u);
    BOOST_TEST(fv[0] == makerow(80u, "abc"));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(one_column_several_rows, RowType, rows_types)
{
    RowType r(1, 42, "abc");
    std::vector<row_view> fv(r.begin(), r.end());
    BOOST_TEST(fv.size() == 2u);
    BOOST_TEST(fv[0] == makerow(42));
    BOOST_TEST(fv[1] == makerow("abc"));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(several_rows_several_columns, RowType, rows_types)
{
    RowType r(2, 80u, "abc", 72u, "cde", 0u, nullptr);
    std::vector<row_view> fv(r.begin(), r.end());
    BOOST_TEST(fv.size() == 3u);
    BOOST_TEST(fv[0] == makerow(80u, "abc"));
    BOOST_TEST(fv[1] == makerow(72u, "cde"));
    BOOST_TEST(fv[2] == makerow(0u, nullptr));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE_TEMPLATE(prefix_increment, RowType, rows_types)
{
    RowType r(2, 80u, "abc", 72u, "cde");
    auto it = r.begin();
    BOOST_TEST(*it == makerow(80u, "abc"));

    // Increment to a dereferenceable state
    BOOST_TEST(&++it == &it);
    BOOST_TEST(*it == makerow(72u, "cde"));

    // Increment to one-past-the-end
    BOOST_TEST(&++it == &it);
    BOOST_TEST(it == r.end());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(postfix_increment, RowType, rows_types)
{
    RowType r(2, 80u, "abc", 72u, "cde");
    auto it = r.begin();
    BOOST_TEST(*it == makerow(80u, "abc"));

    // Increment to a dereferenceable state
    auto itcopy = it++;
    BOOST_TEST(*itcopy == makerow(80u, "abc"));
    BOOST_TEST(*it == makerow(72u, "cde"));

    // Increment to one-past-the-end
    itcopy = it++;
    BOOST_TEST(*itcopy == makerow(72u, "cde"));
    BOOST_TEST(it == r.end());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(prefix_decrement, RowType, rows_types)
{
    RowType r(2, 80u, "abc", 72u, "cde");
    auto it = r.end();

    // Decrement to a dereferenceable state
    BOOST_TEST(&--it == &it);
    BOOST_TEST(*it == makerow(72u, "cde"));

    // Decrement it again
    BOOST_TEST(&--it == &it);
    BOOST_TEST(*it == makerow(80u, "abc"));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(postfix_decrement, RowType, rows_types)
{
    RowType r(2, 80u, "abc", 72u, "cde");
    auto it = r.end();

    // Decrement to a dereferenceable state
    auto itcopy = it--;
    BOOST_TEST(itcopy == r.end());
    BOOST_TEST(*it == makerow(72u, "cde"));

    // Decrement it again
    itcopy = it--;
    BOOST_TEST(*itcopy == makerow(72u, "cde"));
    BOOST_TEST(*it == makerow(80u, "abc"));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(plus_equals, RowType, rows_types)
{
    RowType r(2, 80u, "abc", 72u, "cde", 90u, "fff", 0u, nullptr);
    auto it = r.begin();

    // Increment to a dereferenceable state
    it += 3;
    BOOST_TEST(*it == makerow(0u, nullptr));

    // Increment to one-past-the-end
    it += 1;
    BOOST_TEST(it == r.end());

    // Increment by a negative number
    it += (-2);
    BOOST_TEST(*it == makerow(90u, "fff"));

    // Increment by zero (noop)
    it += 0;
    BOOST_TEST(*it == makerow(90u, "fff"));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(plus_equals_empty, RowType, rows_types)
{
    RowType r;
    auto it = r.begin();
    it += 0;
    BOOST_TEST(it == r.begin());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(minus_equals, RowType, rows_types)
{
    RowType r(2, 80u, "abc", 72u, "cde", 90u, "fff", 0u, nullptr);
    auto it = r.end();

    // Decrement to a dereferenceable state
    it -= 2;
    BOOST_TEST(*it == makerow(90u, "fff"));

    // Decrement to begin
    it -= 2;
    BOOST_TEST(*it == makerow(80u, "abc"));

    // Decrement by a negative number
    it -= (-1);
    BOOST_TEST(*it == makerow(72u, "cde"));

    // Decrement by zero (noop)
    it -= 0;
    BOOST_TEST(*it == makerow(72u, "cde"));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(minus_equals_empty, RowType, rows_types)
{
    RowType r;
    auto it = r.begin();
    it -= 0;
    BOOST_TEST(it == r.begin());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(iterator_plus_ptrdiff, RowType, rows_types)
{
    RowType r(2, 80u, "abc", 72u, "cde", 90u, "fff", 0u, nullptr);
    auto it1 = r.begin();

    // Increment to a dereferenceable state
    auto it2 = it1 + 3;
    BOOST_TEST(*it1 == makerow(80u, "abc"));
    BOOST_TEST(*it2 == makerow(0u, nullptr));

    // Increment to one-past-the-end
    auto it3 = it2 + 1;
    BOOST_TEST(*it2 == makerow(0u, nullptr));
    BOOST_TEST(it3 == r.end());

    // Increment by a negative number
    auto it4 = it3 + (-2);
    BOOST_TEST(it3 == r.end());
    BOOST_TEST(*it4 == makerow(90u, "fff"));

    // Increment by zero (noop)
    auto it5 = it4 + 0;
    BOOST_TEST(*it4 == makerow(90u, "fff"));
    BOOST_TEST(*it5 == makerow(90u, "fff"));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(iterator_plus_ptrdiff_empty, RowType, rows_types)
{
    RowType r;
    BOOST_TEST(r.begin() + 0 == r.begin());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ptrdiff_plus_iterator, RowType, rows_types)
{
    RowType r(2, 80u, "abc", 72u, "cde", 90u, "fff", 0u, nullptr);
    auto it1 = r.begin();

    // Increment to a dereferenceable state
    auto it2 = 3 + it1;
    BOOST_TEST(*it1 == makerow(80u, "abc"));
    BOOST_TEST(*it2 == makerow(0u, nullptr));

    // Increment to one-past-the-end
    auto it3 = 1 + it2;
    BOOST_TEST(*it2 == makerow(0u, nullptr));
    BOOST_TEST(it3 == r.end());

    // Increment by a negative number
    auto it4 = (-2) + it3;
    BOOST_TEST(it3 == r.end());
    BOOST_TEST(*it4 == makerow(90u, "fff"));

    // Increment by zero (noop)
    auto it5 = 0 + it4;
    BOOST_TEST(*it4 == makerow(90u, "fff"));
    BOOST_TEST(*it5 == makerow(90u, "fff"));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ptrdiff_plus_iterator_empty, RowType, rows_types)
{
    RowType r;
    BOOST_TEST(0 + r.begin() == r.begin());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(iterator_minus_ptrdiff, RowType, rows_types)
{
    RowType r(2, 80u, "abc", 72u, "cde", 90u, "fff", 0u, nullptr);
    auto it1 = r.end();

    // Decrement to a dereferenceable state
    auto it2 = it1 - 3;
    BOOST_TEST(it1 == r.end());
    BOOST_TEST(*it2 == makerow(72u, "cde"));

    // Decrement to begin
    auto it3 = it2 - 1;
    BOOST_TEST(*it2 == makerow(72u, "cde"));
    BOOST_TEST(*it3 == makerow(80u, "abc"));

    // Decrement by a negative number
    auto it4 = it3 - (-2);
    BOOST_TEST(*it3 == makerow(80u, "abc"));
    BOOST_TEST(*it4 == makerow(90u, "fff"));

    // Increment by zero (noop)
    auto it5 = it4 - 0;
    BOOST_TEST(*it4 == makerow(90u, "fff"));
    BOOST_TEST(*it5 == makerow(90u, "fff"));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(iterator_minus_ptrdiff_empty, RowType, rows_types)
{
    RowType r;
    BOOST_TEST(r.begin() - 0 == r.begin());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(iterator_minus_iterator, RowType, rows_types)
{
    using It = typename RowType::iterator;

    RowType r(2, 80u, "abc", 72u, "cde", 90u, "fff", 0u, nullptr);
    auto it1 = r.begin();
    auto it2 = r.begin() + 1;
    auto it3 = r.begin() + 2;
    auto it4 = r.begin() + 3;
    auto itend = r.end();

    // Positive
    BOOST_TEST(it2 - it1 == 1);
    BOOST_TEST(it3 - it2 == 1);
    BOOST_TEST(it3 - it1 == 2);
    BOOST_TEST(it4 - it2 == 2);
    BOOST_TEST(itend - it4 == 1);
    BOOST_TEST(itend - it1 == 4);

    // Negative
    BOOST_TEST(it1 - it2 == -1);
    BOOST_TEST(it2 - it3 == -1);
    BOOST_TEST(it1 - it3 == -2);
    BOOST_TEST(it2 - it4 == -2);
    BOOST_TEST(it4 - itend == -1);
    BOOST_TEST(it1 - itend == -4);

    // Zero
    BOOST_TEST(It(it1) - It(it1) == 0);
    BOOST_TEST(It(it2) - It(it2) == 0);
    BOOST_TEST(It(itend) - It(itend) == 0);

    // Self substract
    BOOST_TEST(it1 - it1 == 0);
    BOOST_TEST(it2 - it2 == 0);
    BOOST_TEST(itend - itend == 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(iterator_minus_iterator_empty, RowType, rows_types)
{
    RowType r;
    BOOST_TEST(r.begin() - r.begin() == 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(square_brackets, RowType, rows_types)
{
    RowType r(2, 80u, "abc", 72u, "cde", 90u, "fff", 0u, nullptr);
    auto it = r.begin() + 1;

    BOOST_TEST(it[-1] == makerow(80u, "abc"));
    BOOST_TEST(it[0] == makerow(72u, "cde"));
    BOOST_TEST(it[1] == makerow(90u, "fff"));
    BOOST_TEST(it[2] == makerow(0u, nullptr));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(arrow, RowType, rows_types)
{
    RowType r(2, 80u, "abc", 72u, "cde", 90u, "fff", 0u, nullptr);
    auto it = r.begin() + 1;

    BOOST_TEST(it->size() == 2u);
    BOOST_TEST(it->front() == field_view(72u));
    BOOST_TEST(it->back() == field_view("cde"));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(operator_equals, RowType, rows_types)
{
    using It = typename RowType::iterator;
    RowType empty;
    RowType nonempty(2, 80u, "abc", 72u, "cde");

    struct
    {
        const char* name;
        It it1;
        It it2;
        bool is_equal;
    } test_cases[] = {
        {"value_initialized",               It(),                 It(),                 true },
        {"empty_begin_empty_begin",         empty.begin(),        empty.begin(),        true },
        {"empty_begin_empty_end",           empty.begin(),        empty.end(),          true },
        {"empty_end_empty_end",             empty.end(),          empty.end(),          true },
        {"nonempty_begin_nonempty_begin",   nonempty.begin(),     nonempty.begin(),     true },
        {"nonempty_begin_nonempty_middle",  nonempty.begin(),     nonempty.begin() + 1, false},
        {"nonempty_middle_nonempty_middle", nonempty.begin() + 1, nonempty.begin() + 1, true },
        {"nonempty_begin_nonempty_end",     nonempty.begin(),     nonempty.end(),       false},
        {"nonempty_end_nonempty_end",       nonempty.end(),       nonempty.end(),       true },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            if (tc.is_equal)
            {
                BOOST_TEST(tc.it1 == tc.it2);
                BOOST_TEST(tc.it2 == tc.it1);
                BOOST_TEST(!(tc.it1 != tc.it2));
                BOOST_TEST(!(tc.it2 != tc.it1));
            }
            else
            {
                BOOST_TEST(!(tc.it1 == tc.it2));
                BOOST_TEST(!(tc.it2 == tc.it1));
                BOOST_TEST(tc.it1 != tc.it2);
                BOOST_TEST(tc.it2 != tc.it1);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(operator_lt_lte_gt_gte, RowType, rows_types)
{
    RowType r(2, 80u, "abc", 72u, "cde", 90u, "fff", 0u, nullptr);
    auto it1 = r.begin();
    auto it2 = r.begin() + 1;
    auto it2copy = it2;
    auto it3 = r.begin() + 2;
    auto it4 = r.begin() + 3;
    auto itend = r.end();
    auto itendcopy = itend;

    // it1-it2
    BOOST_TEST(it1 < it2);
    BOOST_TEST(it1 <= it2);
    BOOST_TEST(!(it1 > it2));
    BOOST_TEST(!(it1 >= it2));

    // it2 with itself
    BOOST_TEST(!(it2 < it2copy));
    BOOST_TEST(it2 <= it2copy);
    BOOST_TEST(!(it2 > it2copy));
    BOOST_TEST(it2 >= it2copy);

    // it3-it1
    BOOST_TEST(!(it3 < it1));
    BOOST_TEST(!(it3 <= it1));
    BOOST_TEST(it3 > it1);
    BOOST_TEST(it3 >= it1);

    // it4-itend
    BOOST_TEST(it4 < itend);
    BOOST_TEST(it4 <= itend);
    BOOST_TEST(!(it4 > itend));
    BOOST_TEST(!(it4 >= itend));

    // itend with itself
    BOOST_TEST(!(itend < itendcopy));
    BOOST_TEST(itend <= itendcopy);
    BOOST_TEST(!(itend > itendcopy));
    BOOST_TEST(itend >= itendcopy);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(operator_lt_lte_gt_gte_empty, RowType, rows_types)
{
    RowType r;
    auto it1 = r.begin();
    auto it2 = r.end();

    BOOST_TEST(!(it1 < it2));
    BOOST_TEST(it1 <= it2);
    BOOST_TEST(!(it1 > it2));
    BOOST_TEST(it1 >= it2);
}

BOOST_AUTO_TEST_SUITE_END()

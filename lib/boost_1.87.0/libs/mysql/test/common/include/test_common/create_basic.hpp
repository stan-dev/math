//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_CREATE_BASIC_HPP
#define BOOST_MYSQL_TEST_COMMON_INCLUDE_TEST_COMMON_CREATE_BASIC_HPP

#include <boost/mysql/blob_view.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/row.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/rows.hpp>
#include <boost/mysql/rows_view.hpp>

#include <boost/mysql/detail/access.hpp>

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <vector>

namespace boost {
namespace mysql {
namespace test {

template <class... Types>
BOOST_CXX14_CONSTEXPR std::array<field_view, sizeof...(Types)> make_fv_arr(Types&&... args)
{
    return std::array<field_view, sizeof...(Types)>{{field_view(std::forward<Types>(args))...}};
}

template <class... Types>
std::vector<field_view> make_fv_vector(Types&&... args)
{
    return std::vector<field_view>{field_view(std::forward<Types>(args))...};
}

inline row_view makerowv(const field_view* f, std::size_t size) noexcept
{
    return detail::access::construct<row_view>(f, size);
}

template <class... Types>
row makerow(Types&&... args)
{
    auto fields = make_fv_arr(std::forward<Types>(args)...);
    return row(makerowv(fields.data(), fields.size()));
}

inline rows_view makerowsv(const field_view* fields, std::size_t num_fields, std::size_t num_columns) noexcept
{
    return detail::access::construct<rows_view>(fields, num_fields, num_columns);
}

template <class... Types>
rows makerows(std::size_t num_columns, Types&&... args)
{
    auto fields = make_fv_arr(std::forward<Types>(args)...);
    return rows(makerowsv(fields.data(), fields.size(), num_columns));
}

constexpr time maket(int hours, int mins, int secs, int micros = 0)
{
    return std::chrono::hours(hours) + std::chrono::minutes(mins) + std::chrono::seconds(secs) +
           std::chrono::microseconds(micros);
}

template <std::size_t N>
constexpr string_view makesv(const char (&value)[N])
{
    static_assert(N >= 1, "Expected a C-array literal");
    return string_view(value, N - 1);  // discard null terminator
}

template <std::size_t N>
inline string_view makesv(const std::uint8_t (&value)[N])
{
    return string_view(reinterpret_cast<const char*>(value), N);
}

inline string_view makesv(const std::uint8_t* value, std::size_t size)
{
    return string_view(reinterpret_cast<const char*>(value), size);
}

template <std::size_t N>
blob_view makebv(const char (&value)[N])
{
    static_assert(N >= 1, "Expected a C-array literal");
    return blob_view(reinterpret_cast<const unsigned char*>(value),
                     N - 1);  // discard null terminator
}

template <std::size_t N>
blob makeb(const char (&value)[N])
{
    auto bv = makebv(value);
    return blob(bv.begin(), bv.end());
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif /* TEST_TEST_COMMON_HPP_ */

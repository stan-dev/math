//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_COMMON_TEST_COMMON_HPP
#define BOOST_MYSQL_TEST_COMMON_TEST_COMMON_HPP

#include <boost/mysql/blob_view.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/row.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/rows.hpp>
#include <boost/mysql/rows_view.hpp>

#include <boost/mysql/detail/auxiliar/access_fwd.hpp>
#include <boost/mysql/detail/auxiliar/make_string_view.hpp>
#include <boost/mysql/detail/auxiliar/string_view_offset.hpp>
#include <boost/mysql/detail/protocol/constants.hpp>

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

inline field_view make_svoff_fv(std::size_t offset, std::size_t size, bool is_blob)
{
    return detail::field_view_access::construct(detail::string_view_offset(offset, size), is_blob);
}

inline row_view makerowv(const field_view* f, std::size_t size) noexcept
{
    return detail::row_view_access::construct(f, size);
}

template <class... Types>
row makerow(Types&&... args)
{
    auto fields = make_fv_arr(std::forward<Types>(args)...);
    return row(makerowv(fields.data(), fields.size()));
}

inline rows_view makerowsv(const field_view* fields, std::size_t num_fields, std::size_t num_columns) noexcept
{
    return detail::rows_view_access::construct(fields, num_fields, num_columns);
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
    return detail::make_string_view(value);
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
inline blob_view makebv(const char (&value)[N])
{
    static_assert(N >= 1, "Expected a C-array literal");
    return blob_view(reinterpret_cast<const unsigned char*>(value),
                     N - 1);  // discard null terminator
}

template <std::size_t N>
detail::string_fixed<N> makesfixed(const char (&value)[N + 1])
{
    static_assert(N >= 1, "Expected a C-array literal");
    detail::string_fixed<N> res;
    std::memcpy(res.data(), value, N);
    return res;
}

inline void validate_string_contains(std::string value, const std::vector<std::string>& to_check)
{
    std::transform(value.begin(), value.end(), value.begin(), &tolower);
    for (const auto& elm : to_check)
    {
        BOOST_TEST(
            value.find(elm) != std::string::npos,
            "Substring '" << elm << "' not found in '" << value << "'"
        );
    }
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif /* TEST_TEST_COMMON_HPP_ */

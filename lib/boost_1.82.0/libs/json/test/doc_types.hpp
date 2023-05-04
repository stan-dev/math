//
// Copyright (c) 2022 Dmitry Arkhipov (grisumbras@yandex.ru)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//


#ifndef BOOST_JSON_DOC_DOC_TYPES_HPP
#define BOOST_JSON_DOC_DOC_TYPES_HPP

#if defined(__clang__)
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmismatched-tags"
#endif

#include <array>

//[snippet_conv_spec_trait1
namespace user_ns
{

class ip_address
{
public:
    ip_address(
        unsigned char oct1,
        unsigned char oct2,
        unsigned char oct3,
        unsigned char oct4 );

    const unsigned char*
    begin() const;

    const unsigned char*
    end() const;

private:
    std::array<unsigned char, 4> octets_;
};

template< std::size_t N >
unsigned char
get(const ip_address& addr);

} // namespace user_ns

namespace std
{

template<>
struct tuple_size< user_ns::ip_address >
    : std::integral_constant<std::size_t, 4>
{ };

template< std::size_t N >
struct tuple_element< N, user_ns::ip_address >
{
    using type = unsigned char;
};

} // namespace std
//]

namespace user_ns
{

inline
ip_address::ip_address(
    unsigned char oct1,
    unsigned char oct2,
    unsigned char oct3,
    unsigned char oct4 )
    : octets_{{oct1, oct2, oct3, oct4}}
{ }

inline
const unsigned char*
ip_address::begin() const
{
    return octets_.data();
}

inline
const unsigned char*
ip_address::end() const
{
    return begin() + octets_.size();
}

template< std::size_t N >
unsigned char
get(ip_address const& addr )
{
    static_assert( N >= 0 && N <= 4, "incorrect N" );
    return *(addr.begin() + N);
}

} // namespace user_ns

#endif // BOOST_JSON_DOC_DOC_TYPES_HPP

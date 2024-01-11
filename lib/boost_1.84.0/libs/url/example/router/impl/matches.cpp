//
// Copyright (c) 2023 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#include "../matches.hpp"

namespace boost {
namespace urls {

auto
matches_base::
at( size_type pos ) const
    -> const_reference
{
    if (pos < size())
    {
        return matches()[pos];
    }
    boost::throw_exception(
        std::out_of_range(""));
}

auto
matches_base::
operator[]( size_type pos ) const
    -> const_reference
{
    BOOST_ASSERT(pos < size());
    return matches()[pos];
}

auto
matches_base::
at( core::string_view id ) const
    -> const_reference
{
    for (std::size_t i = 0; i < size(); ++i)
    {
        if (ids()[i] == id)
            return matches()[i];
    }
    boost::throw_exception(
        std::out_of_range(""));
}

auto
matches_base::
operator[]( core::string_view id ) const
    -> const_reference
{
    return at(id);
}

auto
matches_base::
find( core::string_view id ) const
    -> const_iterator
{
    for (std::size_t i = 0; i < size(); ++i)
    {
        if (ids()[i] == id)
            return matches() + i;
    }
    return matches() + size();
}

auto
matches_base::
begin() const
    -> const_iterator
{
    return &matches()[0];
}

auto
matches_base::
end() const
    -> const_iterator
{
    return &matches()[size()];
}

auto
matches_base::
empty() const noexcept
    -> bool
{
    return size() == 0;
}

} // urls
} // boost


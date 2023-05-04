//
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#ifndef BOOST_URL_MATCHES_HPP
#define BOOST_URL_MATCHES_HPP

#include <boost/url/detail/config.hpp>

namespace boost {
namespace urls {

// Base route match results
class matches_base
{
public:
    using iterator = string_view*;
    using const_iterator = string_view const*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference = string_view&;
    using const_reference = string_view const&;
    using pointer = string_view*;
    using const_pointer = string_view const*;

    matches_base() = default;

    virtual ~matches_base() = default;

    virtual
    string_view const*
    matches() const = 0;

    virtual
    string_view const*
    ids() const = 0;

    virtual
    string_view*
    matches() = 0;

    virtual
    string_view*
    ids() = 0;

    virtual
    std::size_t
    size() const = 0;

    virtual
    void
    resize(std::size_t) = 0;

    const_reference
    at( size_type pos ) const;

    const_reference
    at( string_view id ) const;

    const_reference
    operator[]( size_type pos ) const;

    const_reference
    operator[]( string_view id ) const;

    const_iterator
    find( string_view id ) const;

    const_iterator
    begin() const;

    const_iterator
    end() const;

    bool
    empty() const noexcept;
};

/// A range type with the match results
template <std::size_t N = 20>
class matches_storage
    : public matches_base
{
    string_view matches_storage_[N];
    string_view ids_storage_[N];
    std::size_t size_;

    matches_storage(
        string_view matches[N],
        string_view ids[N],
        std::size_t n)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            matches_storage_[i] = matches[i];
            ids_storage_[i] = ids[i];
        }
    }

    virtual
    string_view*
    matches() override
    {
        return matches_storage_;
    }

    virtual
    string_view*
    ids() override
    {
        return ids_storage_;
    }

public:
    matches_storage() = default;

    virtual
    string_view const*
    matches() const override
    {
        return matches_storage_;
    }

    virtual
    string_view const*
    ids() const override
    {
        return ids_storage_;
    }

    virtual
    std::size_t
    size() const override
    {
        return size_;
    }

    virtual
    void
    resize(std::size_t n) override
    {
        size_ = n;
    }
};

/// Default type for storing route match results
using matches = matches_storage<20>;

} // urls
} // boost

#endif


//
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#ifndef BOOST_URL_ROUTER_HPP
#define BOOST_URL_ROUTER_HPP

#include <boost/url/detail/config.hpp>
#include <boost/url/parse_path.hpp>
#include "detail/router.hpp"
#include "matches.hpp"

namespace boost {
namespace urls {

/** A URL router.

    This container matches static and dynamic
    URL requests to an object which represents
    how the it should be handled. These
    values are usually callback functions.

    @tparam T type of resource associated with
    each path template

    @tparam N maximum number of replacement fields
    in a path template

    @par Exception Safety

    @li Functions marked `noexcept` provide the
    no-throw guarantee, otherwise:

    @li Functions which throw offer the strong
    exception safety guarantee.

    @see
        @ref parse_absolute_uri,
        @ref parse_relative_ref,
        @ref parse_uri,
        @ref parse_uri_reference,
        @ref resolve.
*/
template <class T>
class router
    : private detail::router_base
{
public:
    /// Constructor
    router() = default;

    /** Route the specified URL path to a resource

        @param path A url path with dynamic segments
        @param resource A resource the path corresponds to

        @see
            https://fmt.dev/latest/syntax.html
     */
    template <class U>
    void
    insert(core::string_view pattern, U&& v);

    /** Match URL path to corresponding resource

        @param request Request path
        @return The match results
     */
    T const*
    find(segments_encoded_view path, matches_base& m) const noexcept;

#ifdef BOOST_URL_DOCS
    /// @copydoc find
    T const*
    find(segments_encoded_view path, matches& m) const noexcept;
#endif
};

} // urls
} // boost

#include "impl/router.hpp"

#endif


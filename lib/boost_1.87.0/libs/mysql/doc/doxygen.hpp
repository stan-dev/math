//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_DOC_DOXYGEN_HPP
#define BOOST_MYSQL_DOC_DOXYGEN_HPP

/**
 * \brief brief
 * \details
 * description
 *
 * \par Preconditions
 * Only for functions, if they have a precondition (i.e. assert(xxx))
 *
 * \par Exception safety
 * No-throw guarantee.
 * Only for functions. Don't include it in network operations for now.
 * Strong guarantee.
 * Basic guarantee.
 * No-throw guarantee.
 * \throws {exception} {condition}.
 *
 * \par Object lifetimes
 * For async stuff or if the function returns views
 *
 * \par Complexity
 * Include it only for functions in containers or where relevant.
 *      Linear in xxxx.
 *
 * (Include only for async ops)
 * \par Handler signature
 * The handler signature for this operation is `void(boost::mysql::error_code)`
 *
 * \par Per-operation cancellation
 * This operation supports per-operation cancellation. <Describe effects>
 * The following `asio::cancellation_type_t` values are supported:
 *
 *   - `asio::cancellation_type_t::terminal`
 *   - `asio::cancellation_type_t::partial`
 *   - `asio::cancellation_type_t::total`
 *
 * Specify this where it adds any value.
 * \par Thread safety
 * Distinct objects: safe. \n
 * Shared objects: unsafe. \n
 *
 * Specify this for new async operations. Include the error codes that we may return.
 * \par Errors
 * \li \ref client_errc::X when condition Y
 */

#endif

//
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#include <boost/url/detail/except.hpp>
#include <boost/url/decode_view.hpp>
#include <boost/url/grammar/unsigned_rule.hpp>
#include <boost/mp11/algorithm.hpp>

namespace boost {
namespace urls {

template <class T>
template <class U>
void
router<T>::
insert(core::string_view pattern, U&& v)
{
    BOOST_STATIC_ASSERT(
        std::is_same<T, U>::value        ||
        std::is_convertible<U, T>::value ||
        std::is_base_of<T, U>::value);
    using U_ = typename std::decay<
        typename std::conditional<
            std::is_base_of<T, U>::value, U, T
            >::type>::type;
    struct impl : any_resource
    {
        U_ u;

        explicit
            impl(U&& u_)
            : u(std::forward<U>(u_))
        {
        }

        void const*
        get() const noexcept override
        {
            return static_cast<T const*>(&u);
        }
    };
    any_resource const* p = new impl(
        std::forward<U>(v));
    insert_impl( pattern, p );
}

template <class T>
T const*
router<T>::
find(segments_encoded_view path, matches_base& m) const noexcept
{
    core::string_view* matches_it = m.matches();
    core::string_view* ids_it = m.ids();
    any_resource const* p = find_impl(
        path, matches_it, ids_it );
    if (p)
    {
        BOOST_ASSERT(matches_it >= m.matches());
        m.resize(static_cast<std::size_t>(
            matches_it - m.matches()));
        return reinterpret_cast<
            T const*>(p->get());
    }
    m.resize(0);
    return nullptr;
}

} // urls
} // boost

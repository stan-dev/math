/*=============================================================================
    Copyright (c) 2022 niXman (github dot nixman at pm.me)
    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/

#include <boost/fusion/support/pair.hpp>

struct noncopyable_type {
    noncopyable_type(const noncopyable_type &) = delete;
    noncopyable_type& operator=(const noncopyable_type &) = delete;
#if !defined(BOOST_NO_CXX11_RVALUE_REFERENCES)
    noncopyable_type(noncopyable_type &&) = default;
    noncopyable_type& operaotr=(noncopyable_type &&) = default;
#endif
};

int main() {
    using namespace boost::fusion;

    pair<int, noncopyable_type> val = make_pair<int>(noncopyable_type{});

    return 0;
}

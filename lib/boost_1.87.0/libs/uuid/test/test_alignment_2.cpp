// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/uuid.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config.hpp>
#include <boost/config/workaround.hpp>
#include <cstddef>

using namespace boost::uuids;

int main()
{
#if BOOST_WORKAROUND(BOOST_GCC, < 40900)

    BOOST_TEST_LE( alignof(uuid), alignof(max_align_t) );

#else

    BOOST_TEST_LE( alignof(uuid), alignof(std::max_align_t) );

#endif

#if defined(__STDCPP_DEFAULT_NEW_ALIGNMENT__)

    BOOST_TEST_LE( alignof(uuid), __STDCPP_DEFAULT_NEW_ALIGNMENT__ );

#endif

    return boost::report_errors();
}

//  Copyright (c) 2020 Andrey Semashev
//
//  Distributed under the Boost Software License, Version 1.0.
//  See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include <boost/atomic/ipc_atomic_ref.hpp>

#include <boost/config.hpp>
#include <boost/cstdint.hpp>

#include "api_test_helpers.hpp"

int main(int, char *[])
{
    test_lock_free_integral_api< ipc_atomic_ref_wrapper, boost::uint8_t >();
    test_lock_free_integral_api< ipc_atomic_ref_wrapper, boost::int8_t >();

    test_lock_free_integral_api< ipc_atomic_ref_wrapper, boost::uint16_t >();
    test_lock_free_integral_api< ipc_atomic_ref_wrapper, boost::int16_t >();

    test_lock_free_integral_api< ipc_atomic_ref_wrapper, boost::uint32_t >();
    test_lock_free_integral_api< ipc_atomic_ref_wrapper, boost::int32_t >();

    test_lock_free_integral_api< ipc_atomic_ref_wrapper, boost::uint64_t >();
    test_lock_free_integral_api< ipc_atomic_ref_wrapper, boost::int64_t >();

#if defined(BOOST_HAS_INT128) && !defined(BOOST_ATOMIC_TESTS_NO_INT128)
    test_lock_free_integral_api< ipc_atomic_ref_wrapper, boost::int128_type >();
    test_lock_free_integral_api< ipc_atomic_ref_wrapper, boost::uint128_type >();
#endif

#if !defined(BOOST_ATOMIC_NO_FLOATING_POINT)
    test_lock_free_floating_point_api< ipc_atomic_ref_wrapper, float >();
    test_lock_free_floating_point_api< ipc_atomic_ref_wrapper, double >();
    test_lock_free_floating_point_api< ipc_atomic_ref_wrapper, long double >();
#if defined(BOOST_HAS_FLOAT128) && !defined(BOOST_ATOMIC_TESTS_NO_FLOAT128)
    test_lock_free_floating_point_api< ipc_atomic_ref_wrapper, boost::float128_type >();
#endif
#endif

    test_lock_free_pointer_api< ipc_atomic_ref_wrapper, int >();

    test_lock_free_enum_api< ipc_atomic_ref_wrapper >();

    return boost::report_errors();
}

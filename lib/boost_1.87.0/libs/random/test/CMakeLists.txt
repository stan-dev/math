# Copyright 2018, 2019 Peter Dimov
# Copyright 2023 Matt Borland
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

include(BoostTestJamfile OPTIONAL RESULT_VARIABLE HAVE_BOOST_TEST)

if(HAVE_BOOST_TEST)

    boost_test_jamfile(FILE Jamfile.v2 LINK_LIBRARIES
            Boost::random
            Boost::array
            Boost::assert
            Boost::config
            Boost::core
            Boost::dynamic_bitset
            Boost::exception
            Boost::integer
            Boost::io
            Boost::lexical_cast
            Boost::math
            Boost::multiprecision
            Boost::numeric_conversion
            Boost::range
            Boost::static_assert
            Boost::system
            Boost::throw_exception
            Boost::type_traits
            Boost::utility)

endif()

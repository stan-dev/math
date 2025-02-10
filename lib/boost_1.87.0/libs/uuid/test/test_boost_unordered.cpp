// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
#include <boost/core/lightweight_test.hpp>

int main()
{
    using namespace boost::uuids;

    std::size_t const N = 256;

    {
        boost::unordered_flat_set<uuid> set;

        random_generator gen;

        for( std::size_t i = 0; i < N; ++i )
        {
            set.insert( gen() );
        }

        BOOST_TEST_EQ( set.size(), N );
    }

    {
        boost::unordered_flat_set<uuid> set;

        for( std::size_t i = 0; i < N; ++i )
        {
            uuid id = { { static_cast<std::uint8_t>( i ), 0 } };
            set.insert( id );
        }

        BOOST_TEST_EQ( set.size(), N );
    }

    return boost::report_errors();
}

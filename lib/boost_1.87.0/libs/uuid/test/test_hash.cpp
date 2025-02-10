// Copyright (c) 2018 James E. King III
// Copyright (c) 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/container_hash/hash.hpp>
#include <boost/core/lightweight_test.hpp>

int main()
{
    using namespace boost::uuids;

    random_generator gen;

    for( int i = 0; i < 256; ++i )
    {
        uuid u = gen();

        std::size_t v = hash_value( u );

        BOOST_TEST_EQ( boost::hash<uuid>()( u ), v );
        BOOST_TEST_EQ( std::hash<uuid>()( u ), v );
    }

    return boost::report_errors();
}

// Copyright 2010 Andy Tompkins
// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/uuid/namespaces.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/string_generator.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::uuids;

int main()
{
    BOOST_TEST_EQ( ns::dns(), string_generator()("6ba7b810-9dad-11d1-80b4-00c04fd430c8"));
    BOOST_TEST_EQ( ns::url(), string_generator()("6ba7b811-9dad-11d1-80b4-00c04fd430c8"));
    BOOST_TEST_EQ( ns::oid(), string_generator()("6ba7b812-9dad-11d1-80b4-00c04fd430c8"));
    BOOST_TEST_EQ( ns::x500dn(), string_generator()("6ba7b814-9dad-11d1-80b4-00c04fd430c8"));

    return boost::report_errors();
}

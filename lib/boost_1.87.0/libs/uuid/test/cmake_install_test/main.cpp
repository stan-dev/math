// Copyright 2024 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt)

#include <boost/uuid.hpp>

int main()
{
    boost::uuids::time_generator_v1 gen;

    boost::uuids::uuid u1 = gen();
    boost::uuids::uuid u2 = gen();

    return u1 == u2;
}

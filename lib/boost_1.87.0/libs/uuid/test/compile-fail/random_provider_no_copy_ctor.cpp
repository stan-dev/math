// (c) Copyright Andrey Semashev 2018
// (c) Copyright Peter Dimov 2024

// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// https://www.boost.org/LICENSE_1_0.txt)

// The test verifies that random_generator is not copy constructible

#include <boost/uuid/detail/random_provider.hpp>

int main()
{
    boost::uuids::detail::random_provider prov1;
    boost::uuids::detail::random_provider prov2(prov1);

    return 1;
}

// Copyright 2020 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/smart_ptr.hpp>

int main()
{
    boost::shared_ptr<int> p1( new int );
    boost::shared_ptr<int[]> p2( new int[1] );

    boost::make_shared<int>();
    boost::make_shared<int[]>( 1 );
}

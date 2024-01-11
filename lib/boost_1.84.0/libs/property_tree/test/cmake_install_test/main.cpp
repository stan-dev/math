// Copyright 2017, 2021 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// See library home page at http://www.boost.org/libs/system

#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

int main()
{
    pt::ptree tree;

    tree.put( "source.file", __FILE__ );
    tree.put( "source.line", __LINE__ );
}

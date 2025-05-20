// Copyright 2023 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

int main()
{
    pt::ptree tree;

    tree.put( "source.file", __FILE__ );
    tree.put( "source.line", __LINE__ );
}

// Boost.Geometry (aka GGL, Generic Geometry Library)
// QuickBook Example

// Copyright (c) 2011-2024 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//[ring_type
//`Shows how to use the ring_type metafunction, as well as interior_type

#include <iostream>
#include <typeinfo>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

int main()
{
    using point = boost::geometry::model::d2::point_xy<double>;
    using polygon = boost::geometry::model::polygon<point>;

    using ring_type = boost::geometry::ring_type<polygon>::type;
    using int_type = boost::geometry::interior_type<polygon>::type;

    std::cout << typeid(ring_type).name() << std::endl;
    std::cout << typeid(int_type).name() << std::endl;

    // So int_type defines a collection of rings,
    // which is a Boost.Range compatible range
    // The type of an element of the collection is the very same ring type again.
    // We show that.
    using int_ring_type = boost::range_value<int_type>::type;

    std::cout
        << std::boolalpha
        << boost::is_same<ring_type, int_ring_type>::value
        << std::endl;

    return 0;
}

//]

//[ring_type_output
/*`
Output (using gcc):
[pre
          N5boost8geometry5model4ringINS1_2d28point_xyIdNS0_2cs9cartesianEEELb1ELb1ESt6vectorSaEE
St6vectorIN5boost8geometry5model4ringINS2_2d28point_xyIdNS1_2cs9cartesianEEELb1ELb1ES_SaEESaIS9_EE
true
]
*/
//]

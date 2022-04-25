// Boost.Geometry
// Unit Test

// Copyright (c) 2022 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <geometry_test_common.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/helper_geometry.hpp>

class const_point
{
public:
  const_point(double x0, double y0) : x(x0), y(y0) {}

  double GetX() const { return x; }
  double GetY() const { return y; }
private:
  double x{0}, y{0};
};

BOOST_GEOMETRY_REGISTER_POINT_2D_CONST(const_point, double, boost::geometry::cs::cartesian, GetX(), GetY());

using ring_of_const_point = std::vector<const_point>;

// Register a vector of const_pos as a non-const-ring with const points
namespace boost { namespace geometry { namespace traits {
    template<> struct tag<ring_of_const_point> { typedef ring_tag type; };

}}}

template <std::size_t Dim, typename Geometry>
void test_sectionalize_on_const(Geometry const& geometry, std::size_t expected_section_count)
{
    namespace bg = boost::geometry;

    using point_type = typename bg::point_type<Geometry>::type;
    using writable_point_type = typename bg::helper_geometry<point_type>::type;
    using box_type = bg::model::box<writable_point_type>;
    using section_type = bg::sections<box_type, 1>;

    using dim_type = std::integer_sequence<std::size_t, Dim>;

    section_type sections;
    bg::sectionalize<false, dim_type>(geometry, bg::detail::no_rescale_policy(), sections);
    BOOST_CHECK_EQUAL(sections.size(), expected_section_count);
}

int test_main(int, char* [])
{
    ring_of_const_point shape;

    // Construct a shape with 2 monotonic sections in x-dimension and 5 in y-dimension
    shape.emplace_back(0, 0);
    shape.emplace_back(4, 4);
    shape.emplace_back(6, 2);
    shape.emplace_back(8, 4);
    shape.emplace_back(12, 0);
    shape.emplace_back(0, 0);

    // Check if sectionalize (and some other functions) accept this geometry type with const points
    BOOST_CHECK_CLOSE(28, boost::geometry::area(shape), 0.01);
    BOOST_CHECK_CLOSE(28.97, boost::geometry::perimeter(shape), 0.01);
    test_sectionalize_on_const<0>(shape, 2);
    test_sectionalize_on_const<1>(shape, 5);

    return 0;
}


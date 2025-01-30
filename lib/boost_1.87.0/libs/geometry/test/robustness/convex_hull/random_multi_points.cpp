// Boost.Geometry (aka GGL, Generic Geometry Library)
// Robustness Test - convex_hull

// Copyright (c) 2012-2015 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_GEOMETRY_NO_BOOST_TEST

#include <chrono>
#include <sstream>
#include <fstream>

#include <boost/program_options.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/io/svg/svg_mapper.hpp>

#include <common/make_random_generator.hpp>

struct settings_type
{
    bool svg{false};
    bool wkt{false};
    bool verbose{false};
};

namespace bg = boost::geometry;

template <typename Geometry1, typename Geometry2>
void create_svg(std::string const& filename, Geometry1 const& points, Geometry2 const& hull)
{
    typedef typename boost::geometry::point_type<Geometry1>::type point_type;

    boost::geometry::model::box<point_type> box;
    bg::envelope(hull, box);
    bg::buffer(box, box, 1.0);

    std::ofstream svg(filename.c_str());
    boost::geometry::svg_mapper<point_type> mapper(svg, 800, 800);
    mapper.add(box);

    mapper.map(hull, "opacity:0.8;fill:none;stroke:rgb(255,0,255);stroke-width:4;stroke-dasharray:1,7;stroke-linecap:round");
    mapper.map(points, "fill-opacity:0.5;fill:rgb(0,0,255);", 5);
}


template <typename MultiPoint, typename Generator>
inline void make_multi_point(MultiPoint& mp, Generator& generator, int pcount)
{
    typedef typename bg::point_type<MultiPoint>::type point_type;
    typedef typename bg::coordinate_type<MultiPoint>::type coordinate_type;

    for(int i = 0; i < pcount; i++)
    {
        coordinate_type x, y;
        x = generator();
        y = generator();

        point_type p;
        bg::set<0>(p, x);
        bg::set<1>(p, y);

        mp.push_back(p);
    }
}

template <typename MultiPoint, typename Polygon>
bool check_hull(MultiPoint const& mp, Polygon const& poly)
{
    for(typename boost::range_iterator<MultiPoint const>::type it = boost::begin(mp);
        it != boost::end(mp);
        ++it)
    {
        if (! bg::covered_by(*it, poly))
        {
            return false;
        }
    }
    return true;
}


template <typename MultiPoint, typename Generator>
void test_random_multi_points(MultiPoint& result, int& index,
            Generator& generator,
            int pcount, settings_type const& settings)
{
    typedef typename bg::point_type<MultiPoint>::type point_type;

    MultiPoint mp;
    bg::model::polygon<point_type> hull;

    make_multi_point(mp, generator, pcount);
    bg::convex_hull(mp, hull);
    // Check if each point lies in the hull
    bool correct = check_hull(mp, hull);
    if (! correct)
    {
        std::cout << "ERROR! " << std::endl
            << bg::wkt(mp) << std::endl
            << bg::wkt(hull) << std::endl
            << std::endl;
            ;
    }

    if (settings.verbose)
    {
        std::cout << " [" << bg::area(hull) << "]";
    }

    if (settings.svg || ! correct)
    {
        std::ostringstream out;
        out << "random_mp_" << index++ << "_" << pcount << ".svg";
        create_svg(out.str(), mp, hull);
    }
    if (settings.wkt)
    {
        std::cout
            << "input: " << bg::wkt(mp) << std::endl
            << "output: " << bg::wkt(hull) << std::endl
            << std::endl;
            ;
    }
}


template <typename T>
void test_all(int seed, int count, int field_size, int pcount, settings_type const& settings)
{
    auto coordinate_generator = make_int_generator(seed, field_size - 1);

    typedef bg::model::multi_point
        <
            bg::model::d2::point_xy<T>
        > mp;


    auto const t0 = std::chrono::high_resolution_clock::now();
    int index = 0;
    for(int i = 0; i < count; i++)
    {
        mp p;
        test_random_multi_points<mp>(p, index, coordinate_generator, pcount, settings);
    }
    auto const t = std::chrono::high_resolution_clock::now();
    auto const elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t - t0).count();
    std::cout << std::endl
        << "points: " << index
        << " type: " << typeid(T).name()
        << " time: " << elapsed_ms / 1000.0  << std::endl;
}

int main(int argc, char** argv)
{
    try
    {
        namespace po = boost::program_options;
        po::options_description description("=== random_multi_points ===\nAllowed options");

        std::string type = "double";
        int count = 1;
        int seed = -1;
        int pcount = 30;
        int field_size = 10;
        settings_type settings;

        description.add_options()
            ("help", "Help message")
            ("seed", po::value<int>(&seed), "Initialization seed for random generator")
            ("count", po::value<int>(&count), "Number of tests")
            ("number", po::value<int>(&pcount), "Number of points")
            ("size", po::value<int>(&field_size), "Size of the field")
            ("type", po::value<std::string>(&type)->default_value("double"), "Type (int,float,double)")
            ("verbose", po::value<bool>(&settings.verbose), "Verbose")
            ("wkt", po::value<bool>(&settings.wkt), "Create a WKT of the inputs, for all tests")
            ("svg", po::value<bool>(&settings.svg), "Create a SVG for all tests")
        ;

        po::variables_map varmap;
        po::store(po::parse_command_line(argc, argv, description), varmap);
        po::notify(varmap);

        if (varmap.count("help"))
        {
            std::cout << description << std::endl;
            return 1;
        }

        if (type == "float")
        {
            test_all<float>(seed, count, field_size, pcount, settings);
        }
        else if (type == "double")
        {
            test_all<double>(seed, count, field_size, pcount, settings);
        }
        else if (type == "int")
        {
            test_all<int>(seed, count, field_size, pcount, settings);
        }

    }
    catch(std::exception const& e)
    {
        std::cout << "Exception " << e.what() << std::endl;
    }
    catch(...)
    {
        std::cout << "Other exception" << std::endl;
    }

    return 0;
}

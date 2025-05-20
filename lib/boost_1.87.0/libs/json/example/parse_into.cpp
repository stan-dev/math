//
// Copyright (c) 2021 Peter Dimov
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

//
// An example that compares the performance of json::parse and
// json::parse_into, using https://github.com/kostya/benchmarks#json
//
// Typical results:
//
// 1.json: 115075437 bytes
// boost::json::parse: 721 ms
//   x: -5.00335e-30, y: 5.00428e+30, z: 0.499722: 121 ms
// parse_into coordinates: 398 ms
//   x: -5.00335e-30, y: 5.00428e+30, z: 0.499722: 3 ms
// parse_into coordinates2: 326 ms
//   x: -5.00335e-30, y: 5.00428e+30, z: 0.499722: 0 ms
//

#include <boost/json.hpp>
#include <iostream>

#if !defined(BOOST_DESCRIBE_CXX14)

#include <boost/config/pragma_message.hpp>

BOOST_PRAGMA_MESSAGE( "This example requires C++14" )

int main() {}

#else

#include <boost/describe.hpp>

#include <chrono>
#include <fstream>
#include <iterator>
#include <map>
#include <vector>

// An std::map<std::string, std::pair<int, bool>> replacement
// We don't need to store the options

struct options
{
    using mapped_type = std::pair<int, bool>;
    using value_type = std::pair<std::string, mapped_type>;
    using iterator = value_type*;

    std::pair<iterator, bool>
    emplace( value_type const& );

    void
    emplace( std::string const&, mapped_type const& )
    {
    }

    iterator
    begin();

    iterator
    end();

    void
    clear()
    { }
};

struct coordinate
{
    double x{}, y{}, z{};
    std::string name;
    options opts;
};

BOOST_DESCRIBE_STRUCT(coordinate, (), (x, y, z, name, opts))

struct coordinates1
{
    std::vector<coordinate> coordinates;
    std::string info;
};

BOOST_DESCRIBE_STRUCT(coordinates1, (), (coordinates, info))

// std::vector<coordinate> replacement that just
// keeps a running sum

struct accumulator
{
    using value_type = coordinate;
    using iterator = coordinate*;

    std::size_t len = 0;

    double x = 0;
    double y = 0;
    double z = 0;

    void push_back( coordinate const& v )
    {
        x += v.x;
        y += v.y;
        z += v.z;

        ++len;
    }

    iterator
    begin() { return nullptr; }

    iterator
    end() { return nullptr; }

    void
    clear()
    { }
};

struct coordinates2
{
    accumulator coordinates;
    std::string info;
};

BOOST_DESCRIBE_STRUCT(coordinates2, (), (coordinates, info))

using namespace std::chrono_literals;

int main()
{
    // https://github.com/kostya/benchmarks/blob/master/json/generate_json.rb
    std::ifstream is( "/tmp/1.json" );
    std::string json( std::istreambuf_iterator<char>( is ), std::istreambuf_iterator<char>{} );

    std::cout << "1.json: " << json.size() << " bytes\n";

    // https://github.com/kostya/benchmarks/blob/master/json/test_boost_json.cpp
    {
        auto tp1 = std::chrono::steady_clock::now();

        boost::json::value jv = boost::json::parse( json );

        auto tp2 = std::chrono::steady_clock::now();
        std::cout << "boost::json::parse: " << (tp2 - tp1) / 1ms << " ms\n";

        auto x = 0.0, y = 0.0, z = 0.0;
        auto len = 0;

        auto &obj = jv.get_object();

        for( auto& v: obj["coordinates"].get_array() )
        {
            ++len;
            auto& coord = v.get_object();
            x += coord["x"].get_double();
            y += coord["y"].get_double();
            z += coord["z"].get_double();
        }

        x /= len;
        y /= len;
        z /= len;

        auto tp3 = std::chrono::steady_clock::now();
        std::cout << "  x: " << x << ", y: " << y << ", z: " << z << ": " << (tp3 - tp2) / 1ms << " ms\n";
    }

    {
        auto tp1 = std::chrono::steady_clock::now();

        coordinates1 w;

        boost::system::error_code ec;
        boost::json::parse_into( w, json, ec );

        if( ec.failed() )
        {
            std::cout << "Error: " << ec.what() << std::endl;
        }

        auto tp2 = std::chrono::steady_clock::now();
        std::cout << "parse_into coordinates: " << (tp2 - tp1) / 1ms << " ms\n";

        auto x = 0.0, y = 0.0, z = 0.0;
        auto len = 0;

        for( auto const& v: w.coordinates )
        {
            x += v.x;
            y += v.y;
            z += v.z;

            ++len;
        }

        x /= len;
        y /= len;
        z /= len;

        auto tp3 = std::chrono::steady_clock::now();
        std::cout << "  x: " << x << ", y: " << y << ", z: " << z << ": " << (tp3 - tp2) / 1ms << " ms\n";
    }

    {
        auto tp1 = std::chrono::steady_clock::now();

        coordinates2 w;

        boost::system::error_code ec;
        boost::json::parse_into( w, json, ec );

        if( ec.failed() )
        {
            std::cout << "Error: " << ec.what() << std::endl;
        }

        auto tp2 = std::chrono::steady_clock::now();
        std::cout << "parse_into coordinates2: " << (tp2 - tp1) / 1ms << " ms\n";

        double x = w.coordinates.x / w.coordinates.len;
        double y = w.coordinates.y / w.coordinates.len;
        double z = w.coordinates.z / w.coordinates.len;

        auto tp3 = std::chrono::steady_clock::now();
        std::cout << "  x: " << x << ", y: " << y << ", z: " << z << ": " << (tp3 - tp2) / 1ms << " ms\n";
    }
}

#endif

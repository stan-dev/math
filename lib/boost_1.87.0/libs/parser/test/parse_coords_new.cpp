/**
 *   Copyright (C) 2024 Phil Endecott
 *   Copyright (C) 2024 T. Zachary Laine
 *
 *   Distributed under the Boost Software License, Version 1.0. (See
 *   accompanying file LICENSE_1_0.txt or copy at
 *   http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/parser/config.hpp>

#include <boost/core/lightweight_test.hpp>

#if BOOST_PARSER_USE_STD_TUPLE
int main() { return boost::report_errors(); }
#else

#include <boost/parser/parser.hpp>

#include <cassert>
#include <iostream>
#include <optional>
#include <string>


namespace bp = boost::parser;

using namespace boost::hana::literals;

namespace g2d {
    struct Vector
    {
        double x, y;
    };

    std::ostream & operator<<(std::ostream & os, Vector vec)
    {
        os << "vec{" << vec.x << ", " << vec.y << "}";
        return os;
    }
};

bp::rule<struct uint_0_60,               unsigned int> uint_0_60               = "uint_0_60";
bp::rule<struct double_0_60,             double>       double_0_60             = "double_0_60";
bp::rule<struct degrees_decimal_minutes, double>       degrees_decimal_minutes = "degrees_decimal_minutes";
bp::rule<struct degrees_minutes_seconds, double>       degrees_minutes_seconds = "degrees_minutes_seconds";
bp::rule<struct degrees,                 double>       degrees                 = "degrees";
bp::rule<struct latitude,                double>       latitude                = "latitude";
bp::rule<struct longitude,               double>       longitude               = "longitude";
bp::rule<struct signed_latitude,         double>       signed_latitude         = "signed_latitude";
bp::rule<struct signed_longitude,        double>       signed_longitude        = "signed_longitude";
bp::rule<struct latlon,                  g2d::Vector>  latlon                  = "latlon";


const auto degrees_symbol = bp::no_case[                 bp::lit("degrees") | bp::lit("deg") | bp::lit('d') ];
const auto minutes_symbol = bp::no_case[ bp::lit('\'') | bp::lit("minutes") | bp::lit("min") | bp::lit('m') ];
const auto seconds_symbol = bp::no_case[ bp::lit('"')  | bp::lit("seconds") | bp::lit("sec") | bp::lit('s') ];

const auto uint_0_60_def   = bp::uint_   [( [](auto & ctx) { _pass(ctx) = _attr(ctx) < 60U; _val(ctx) = _attr(ctx); } )];
const auto double_0_60_def = bp::double_ [( [](auto & ctx) { _pass(ctx) = _attr(ctx) < 60;  _val(ctx) = _attr(ctx); } )];

const auto decimal_degrees = bp::double_ >> -degrees_symbol;

const auto degrees_decimal_minutes_def = (bp::uint_ >> -degrees_symbol
                                       >> (double_0_60 - '.') >> -minutes_symbol) [( [](auto & ctx) {
  auto d = _attr(ctx)[0_c];
  auto m = _attr(ctx)[1_c];
  _val(ctx) = d + m/60.0;
} )];

const auto degrees_minutes_seconds_def = (bp::uint_ >> -degrees_symbol
                                       >> uint_0_60 >> -minutes_symbol
                                       >> (double_0_60 - '.') >> -seconds_symbol) [( [](auto & ctx) {
  auto d = _attr(ctx)[0_c];
  auto m = _attr(ctx)[1_c];
  auto s = _attr(ctx)[2_c];
  _val(ctx) = d + m/60.0 + s/3600.0;
} )];

const auto degrees_def = degrees_minutes_seconds
                       | degrees_decimal_minutes
                       | decimal_degrees;

const auto northsouth = bp::no_case[ bp::char_("ns") ];
const auto eastwest   = bp::no_case[ bp::char_("ew") ];

const auto latitude_def = (degrees >> northsouth) [( [](auto & ctx) {
  auto d  = _attr(ctx)[0_c];
  auto ns = _attr(ctx)[1_c];
  _pass(ctx) = d <= 90;
  _val(ctx) = ns=='S' || ns=='s' ? -d : d;
} )];

const auto longitude_def = (degrees >> eastwest) [( [](auto & ctx) {
  auto d  = _attr(ctx)[0_c];
  auto ew = _attr(ctx)[1_c];
  _pass(ctx) = d <= 180;
  _val(ctx) = ew=='W' || ew=='w' ? -d : d;
} )];

const auto signed_degrees = bp::double_ >> -degrees_symbol;

const auto signed_latitude_def  = signed_degrees [( [](auto & ctx) { auto d = _attr(ctx); _pass(ctx) =  -90 <= d && d <=  90; _val(ctx) = _attr(ctx); } )];
const auto signed_longitude_def = signed_degrees [( [](auto & ctx) { auto d = _attr(ctx); _pass(ctx) = -180 <= d && d <= 180; _val(ctx) = _attr(ctx); } )];

const auto latlon_def = ((latitude >> longitude) [( [](auto & ctx) { _val(ctx) = g2d::Vector{_attr(ctx)[1_c], _attr(ctx)[0_c]}; } )] )
                      | ((longitude >> latitude) [( [](auto & ctx) { _val(ctx) = g2d::Vector{_attr(ctx)[0_c], _attr(ctx)[1_c]}; } )] )
                      | ((signed_longitude >> signed_latitude)
                                                 [( [](auto & ctx) { _val(ctx) = g2d::Vector{_attr(ctx)[0_c], _attr(ctx)[1_c]}; } )] );

BOOST_PARSER_DEFINE_RULES(uint_0_60);
BOOST_PARSER_DEFINE_RULES(double_0_60);
BOOST_PARSER_DEFINE_RULES(degrees_decimal_minutes);
BOOST_PARSER_DEFINE_RULES(degrees_minutes_seconds);
BOOST_PARSER_DEFINE_RULES(degrees);
BOOST_PARSER_DEFINE_RULES(latitude);
BOOST_PARSER_DEFINE_RULES(longitude);
BOOST_PARSER_DEFINE_RULES(signed_latitude);
BOOST_PARSER_DEFINE_RULES(signed_longitude);
BOOST_PARSER_DEFINE_RULES(latlon);


#if 0
static std::optional<g2d::Vector> try_parse_latlon(std::string_view s)
{
  auto input = latlon >> bp::eoi;
  return parse(s,input, bp::ws|bp::lit(','));
}


int main(int argc, const char* argv[])
{
  assert(argc==2);

  auto opt_coords = try_parse_latlon(argv[1]);

  if (!opt_coords) {
    std::cerr << "parse error\n";
    return 1;
  } else {
    std::cout << opt_coords->x << " " << opt_coords->y << "\n";
    return 0;
  }
}
#endif

int main()
{
    std::vector<std::string> test_coords = {
        "12.34 N, 56.78 E",
        "56.78,-12.34",
        "12d 30' n 45d 15' 7\" w",
        "12 30 45 N, 45 15 7 W",
        "12d 30m 15s S 50d 59m 59s W",
        "50d 0.5m n 50d 59m 59s W"};

    {
        auto result = bp::parse(test_coords[0], latlon, bp::ws | bp::lit(','));
        BOOST_TEST(result);
        BOOST_TEST(std::abs(result->x - 56.78) < 0.001);
        BOOST_TEST(std::abs(result->y - 12.34) < 0.001);
    }

    {
        auto result = bp::parse(test_coords[1], latlon, bp::ws | bp::lit(','));
        BOOST_TEST(result);
        BOOST_TEST(std::abs(result->x - 56.78) < 0.001);
        BOOST_TEST(std::abs(result->y - -12.34) < 0.001);
    }

    {
        auto result = bp::parse(test_coords[2], latlon, bp::ws | bp::lit(','));
        BOOST_TEST(result);
        BOOST_TEST(std::abs(result->x - -45.2519) < 0.001);
        BOOST_TEST(std::abs(result->y - 12.5) < 0.001);
    }

    {
        auto result = bp::parse(test_coords[3], latlon, bp::ws | bp::lit(','));
        BOOST_TEST(result);
        BOOST_TEST(std::abs(result->x - -45.2519) < 0.001);
        BOOST_TEST(std::abs(result->y - 12.5125) < 0.001);
    }

    {
        auto result = bp::parse(test_coords[4], latlon, bp::ws | bp::lit(','));
        BOOST_TEST(result);
        BOOST_TEST(std::abs(result->x - -50.9997) < 0.001);
        BOOST_TEST(std::abs(result->y - -12.5042) < 0.001);
    }

    {
        auto result = bp::parse(test_coords[5], latlon, bp::ws | bp::lit(','));
        BOOST_TEST(result);
        BOOST_TEST(std::abs(result->x - -50.9997) < 0.001);
        BOOST_TEST(std::abs(result->y - 50.0083) < 0.001);
    }

    return boost::report_errors();
}

#endif

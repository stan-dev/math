// Copyright 2019 Przemyslaw Bartosik
// Copyright 2019 Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>
#include <boost/histogram/accumulators/mean.hpp>
#include <boost/histogram/accumulators/ostream.hpp>
#include <boost/histogram/axis/category.hpp>
#include <boost/histogram/axis/integer.hpp>
#include <boost/histogram/axis/option.hpp>
#include <boost/histogram/axis/regular.hpp>
#include <boost/histogram/make_histogram.hpp>
#include <boost/histogram/ostream.hpp>
#include <limits>
#include <sstream>
#include <string>
#include "throw_exception.hpp"
#include "utility_histogram.hpp"

using namespace boost::histogram;

template <class Histogram>
auto str(const Histogram& h, const unsigned width = 0) {
  std::ostringstream os;
  // BEGIN and END make nicer error messages
  os << "BEGIN\n" << std::setw(width) << h << "END";
  return os.str();
}

template <class Tag>
void run_tests() {
  using R = axis::regular<>;
  using R2 =
      axis::regular<double, boost::use_default, axis::null_type, axis::option::none_t>;
  using R3 = axis::regular<double, axis::transform::log>;
  using C = axis::category<std::string>;
  using I = axis::integer<>;

  // regular
  {
    auto h = make(Tag(), R(3, -0.5, 1.0));
    h.at(0) = 1;
    h.at(1) = 10;
    h.at(2) = 5;

    const auto expected =
        "BEGIN\n"
        "histogram(regular(3, -0.5, 1, options=underflow | overflow))\n"
        "                +------------------------------------------------------------+\n"
        "[-inf, -0.5) 0  |                                                            |\n"
        "[-0.5,    0) 1  |======                                                      |\n"
        "[   0,  0.5) 10 |=========================================================== |\n"
        "[ 0.5,    1) 5  |==============================                              |\n"
        "[   1,  inf) 0  |                                                            |\n"
        "                +------------------------------------------------------------+\n"
        "END";

    BOOST_TEST_CSTR_EQ(expected, str(h).c_str());
  }

  // regular, narrow
  {
    auto h = make(Tag(), R2(3, -0.5, 1.0));
    h.at(0) = 1;
    h.at(1) = 10;
    h.at(2) = 2;

    const auto expected = "BEGIN\n"
                          "histogram(regular(3, -0.5, 1, options=none))\n"
                          "               +-----------------------+\n"
                          "[-0.5,   0) 1  |==                     |\n"
                          "[   0, 0.5) 10 |====================== |\n"
                          "[ 0.5,   1) 2  |====                   |\n"
                          "               +-----------------------+\n"
                          "END";

    BOOST_TEST_CSTR_EQ(expected, str(h, 40).c_str());

    // too narrow
    BOOST_TEST_CSTR_EQ("BEGIN\n"
                       "histogram(regular(3, -0.5, 1, options=none))END",
                       str(h, 10).c_str());
  }

  // regular2
  {
    auto h = make(Tag(), R2(3, -0.5, 1.0));
    h.at(0) = 1;
    h.at(1) = -5;
    h.at(2) = 2;

    const auto expected =
        "BEGIN\n"
        "histogram(regular(3, -0.5, 1, options=none))\n"
        "               +-------------------------------------------------------------+\n"
        "[-0.5,   0) 1  |                                           =========         |\n"
        "[   0, 0.5) -5 |===========================================                  |\n"
        "[ 0.5,   1) 2  |                                           ================= |\n"
        "               +-------------------------------------------------------------+\n"
        "END";

    BOOST_TEST_CSTR_EQ(expected, str(h).c_str());
  }

  // regular with log
  {
    auto h = make(Tag(), R3(6, 1e-3, 1e3, "foo"));

    const auto expected =
        "BEGIN\n"
        "histogram(regular_log(6, 0.001, 1000, metadata=\"foo\", options=underflow | "
        "overflow))\n"
        "                 +-----------------------------------------------------------+\n"
        "[    0, 0.001) 0 |                                                           |\n"
        "[0.001,  0.01) 0 |                                                           |\n"
        "[ 0.01,   0.1) 0 |                                                           |\n"
        "[  0.1,     1) 0 |                                                           |\n"
        "[    1,    10) 0 |                                                           |\n"
        "[   10,   100) 0 |                                                           |\n"
        "[  100,  1000) 0 |                                                           |\n"
        "[ 1000,   inf) 0 |                                                           |\n"
        "                 +-----------------------------------------------------------+\n"
        "END";

    BOOST_TEST_CSTR_EQ(expected, str(h).c_str());
  }

  // integer
  {
    auto h = make(Tag(), I(0, 1));
    h.at(0) = -10;
    h.at(1) = 5;

    const auto expected =
        "BEGIN\n"
        "histogram(integer(0, 1, options=underflow | overflow))\n"
        "       +---------------------------------------------------------------------+\n"
        "-1 0   |                                                                     |\n"
        " 0 -10 |=============================================                        |\n"
        " 1 5   |                                             ======================= |\n"
        "       +---------------------------------------------------------------------+\n"
        "END";

    BOOST_TEST_CSTR_EQ(expected, str(h).c_str());
  }

  // catorgy<string>
  {
    auto h = make(Tag(), C({"a", "bb", "ccc", "dddd"}));
    h.at(0) = 1.23;
    h.at(1) = 1;
    h.at(2) = 1.2345789e-3;
    h.at(3) = 1.2345789e-12;
    h.at(4) = std::numeric_limits<double>::quiet_NaN();

    const auto expected =
        "BEGIN\n"
        "histogram(category(\"a\", \"bb\", \"ccc\", \"dddd\", options=overflow))\n"
        "                +------------------------------------------------------------+\n"
        "    a 1.23      |=========================================================== |\n"
        "   bb 1         |================================================            |\n"
        "  ccc 0.001235  |                                                            |\n"
        " dddd 1.235e-12 |                                                            |\n"
        "other nan       |                                                            |\n"
        "                +------------------------------------------------------------+\n"
        "END";

    BOOST_TEST_CSTR_EQ(expected, str(h).c_str());
  }

  // histogram with axis that has no value method
  {
    struct minimal_axis {
      int index(int x) const { return x % 2; }
      int size() const { return 2; }
    };

    auto h = make(Tag(), minimal_axis{});
    h.at(0) = 3;
    h.at(1) = 4;

    const auto expected =
        "BEGIN\n"
        "histogram(<unstreamable>)\n"
        "    +------------------------------------------------------------------------+\n"
        "0 3 |=====================================================                   |\n"
        "1 4 |======================================================================= |\n"
        "    +------------------------------------------------------------------------+\n"
        "END";

    BOOST_TEST_CSTR_EQ(expected, str(h).c_str());
  }

  // fallback for 2D
  {
    auto h = make(Tag(), R(1, -1, 1), R(2, -4, 7));
    h.at(-1, 0) = 1000;
    h.at(-1, -1) = 123;
    h.at(1, 0) = 1.23456789;
    h.at(-1, 2) = std::numeric_limits<double>::quiet_NaN();

    const auto expected =
        "BEGIN\n"
        "histogram(\n"
        "  regular(1, -1, 1, options=underflow | overflow)\n"
        "  regular(2, -4, 7, options=underflow | overflow)\n"
        "  (-1 -1): 123   ( 0 -1): 0     ( 1 -1): 0     (-1  0): 1000 \n"
        "  ( 0  0): 0     ( 1  0): 1.235 (-1  1): 0     ( 0  1): 0    \n"
        "  ( 1  1): 0     (-1  2): nan   ( 0  2): 0     ( 1  2): 0    \n"
        ")END";

    BOOST_TEST_CSTR_EQ(expected, str(h).c_str());
  }

  // fallback for profile
  {
    auto h = make_s(Tag(), profile_storage(), R(1, -1, 1));
    h.at(0) = accumulators::mean<>(10, 100, 1000);

    const auto expected = "BEGIN\n"
                          "histogram(\n"
                          "  regular(1, -1, 1, options=underflow | overflow)\n"
                          "  (-1): mean(0, 0, -0)      ( 0): mean(10, 100, 1000)\n"
                          "  ( 1): mean(0, 0, -0)     \n"
                          ")END";

    BOOST_TEST_CSTR_EQ(expected, str(h).c_str());
  }
}

int main() {
  run_tests<static_tag>();
  run_tests<dynamic_tag>();

  {
    // cannot make empty static histogram
    auto h = histogram<std::vector<axis::regular<>>>();

    const auto expected = "BEGIN\n"
                          "histogram()END";

    BOOST_TEST_CSTR_EQ(expected, str(h).c_str());
  }

  return boost::report_errors();
}

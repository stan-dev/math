//  (C) Copyright Raffi Enficiaud 2019.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
// ***************************************************************************

// Boost.Test
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_parameters.hpp>

bool init_unit_test()
{
  using namespace boost::unit_test;

// Having some problems on AppleClang 10.10 / Xcode 6/7
// Persists on MacOS 11 with AppleClang 13
#if (!defined(BOOST_TEST_DYN_LINK) || (!defined(BOOST_CLANG) || (BOOST_CLANG != 1) || (__clang_major__ >= 8))) && !defined(__APPLE__)
  log_level logLevel = runtime_config::get<log_level>(runtime_config::btrt_log_level);
  std::cout << "Current log level: " << static_cast<int>(logLevel) << std::endl;
  output_format logFormat = runtime_config::get<output_format>(runtime_config::btrt_log_format);
  std::cout << "Current log format: " << static_cast<int>(logFormat) << std::endl;
  report_level reportLevel = runtime_config::get<report_level>(runtime_config::btrt_report_level);
  std::cout << "Current report level: " << static_cast<int>(reportLevel) << std::endl;
#endif
  return true;
}

BOOST_AUTO_TEST_CASE( my_test1 )
{
    BOOST_CHECK( true );
}

int main(int argc, char* argv[])
{
  int retCode = boost::unit_test::unit_test_main( &init_unit_test, argc, argv );
  return retCode;
}

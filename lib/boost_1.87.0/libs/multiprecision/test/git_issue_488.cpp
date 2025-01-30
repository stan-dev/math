///////////////////////////////////////////////////////////////////////////////
//  Copyright 2022 John Maddock. Distributed under the Boost
//  Copyright 2022 Christopher Kormanyos. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>

#include <boost/multiprecision/cpp_int.hpp>

#include "test.hpp"

int main()
{
   using local_uint265_type = boost::multiprecision::uint256_t;

   std::uint16_t bits[16] = { };
   std::uint8_t  result_bits[16] = { };

   std::fill(std::begin(bits), std::end(bits), static_cast<std::uint16_t>(UINT16_C(0x5678)));
   std::fill(std::begin(result_bits), std::end(result_bits), static_cast<std::uint8_t>(bits[0] & 0xFF));

   local_uint265_type u{}, v{};

   const std::array<bool, static_cast<std::size_t>(UINT8_C(2))> msv_values{ true, false };

   const auto flg = std::cout.flags();

   for (const auto msv_first : msv_values)
   {
      static_cast<void>(import_bits(u, std::begin(bits), std::end(bits), 8u, msv_first));
      static_cast<void>(import_bits(v, std::begin(result_bits), std::end(result_bits), 8u, msv_first));

      std::cout << std::hex << std::uppercase << u << std::endl;

      BOOST_CHECK_EQUAL(u, v);
   }

   std::cout.flags(flg);

   return boost::report_errors();
}

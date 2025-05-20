// Copyright 2022 Christian Mazakas
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/unordered/detail/narrow_cast.hpp>

#include <boost/core/lightweight_test.hpp>
#include <boost/cstdint.hpp>

// want to prove that for the wider type, the higher bits of the value
// represenation don't affect the results of the narrowing, which in this case
// is masking out the high bits when comapred to the narrow type

static void signed_integral_narrowing()
{
  // test positive range, fits
  // [0, 127]
  for (boost::int32_t i = 0x00; i < 0x80; ++i) {
    boost::int8_t k = (boost::int8_t)i;
    BOOST_TEST_GE(k, 0);
    BOOST_TEST_EQ(boost::unordered::detail::narrow_cast<boost::int8_t>(i), k);
  }

  // test positive range, doesn't fit
  // [0xff00, 0xff7f]
  for (boost::int32_t i = 0x00; i < 0x80; ++i) {
    boost::int32_t j = i + 0xff00;
    boost::int8_t k = (boost::int8_t)i;
    BOOST_TEST_GE(k, 0);
    BOOST_TEST_EQ(boost::unordered::detail::narrow_cast<boost::int8_t>(j), k);
  }

  // test negative range, fits
  // [-128, -1]
  for (boost::int32_t i = 0x00; i < 0x80; ++i) {
    boost::int32_t j = i + (boost::int32_t)0xffffff80;
    boost::int8_t k = (boost::int8_t)j;
    BOOST_TEST_LT(j, 0);
    BOOST_TEST_LT(k, 0);
    BOOST_TEST_EQ(boost::unordered::detail::narrow_cast<boost::int8_t>(j), k);
  }

  // test negative range, doesn't fit
  for (boost::int32_t i = 0x00; i < 0x80; ++i) {
    boost::int32_t j = i + (boost::int32_t)0x80000000;
    boost::int8_t k = (boost::int8_t)(i);
    BOOST_TEST_LT(j, 0);
    BOOST_TEST_EQ(boost::unordered::detail::narrow_cast<boost::int8_t>(j), k);
  }

  for (boost::int32_t i = 0x00; i < 0x100; ++i) {
    boost::int32_t j = (boost::int32_t)0x80ff0000 + i;
    BOOST_TEST_LT(j, 0);
    BOOST_TEST_EQ(boost::unordered::detail::narrow_cast<boost::int8_t>(j),
      (boost::int8_t)i);
  }

  // test special values
  {
    boost::int32_t x = 0xff;
    BOOST_TEST_EQ(boost::unordered::detail::narrow_cast<boost::int8_t>(x), -1);
  }

  {
    boost::int32_t x = (boost::int32_t)0xffffff00;
    BOOST_TEST_EQ(boost::unordered::detail::narrow_cast<boost::int8_t>(x),
      (boost::int8_t)0x00);
  }

  {
    boost::int32_t x = (boost::int32_t)0xffffff7f;
    BOOST_TEST_EQ(boost::unordered::detail::narrow_cast<boost::int8_t>(x),
      (boost::int8_t)0x7f);
  }

  {
    boost::int32_t x = (boost::int32_t)0xffffffff;
    BOOST_TEST_EQ(boost::unordered::detail::narrow_cast<boost::int8_t>(x),
      (boost::int8_t)-1);
  }
}

static void unsigned_integral_narrowing()
{
  // test range: [0x00, 0xff]
  for (boost::uint32_t i = 0x00; i < 0x100; ++i) {
    BOOST_TEST_EQ(boost::unordered::detail::narrow_cast<boost::uint8_t>(i),
      (boost::uint8_t)(i & 0xff));
  }

  // test range: [0xffffff00, 0xffffffff]
  boost::uint32_t i = 0xffffff00;
  for (; i < 0xffffffff; ++i) {
    BOOST_TEST_EQ(boost::unordered::detail::narrow_cast<boost::uint8_t>(i),
      (boost::uint8_t)(i & 0xff));
  }
  BOOST_TEST_EQ(boost::unordered::detail::narrow_cast<boost::uint8_t>(i),
    (boost::uint8_t)(i & 0xff));
}

int main()
{
  signed_integral_narrowing();
  unsigned_integral_narrowing();

  return boost::report_errors();
}

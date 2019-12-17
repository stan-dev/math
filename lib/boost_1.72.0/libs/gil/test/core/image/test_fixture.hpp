//
// Copyright 2019 Mateusz Loskot <mateusz at loskot dot net>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//
#include <boost/gil.hpp>
#include <boost/assert.hpp>

#include <cstdint>
#include <initializer_list>
#include <limits>
#include <random>
#include <tuple>
#include <type_traits>

namespace boost { namespace gil {

namespace test { namespace fixture {

using image_types = std::tuple
<
    gil::gray8_image_t,
    gil::gray16_image_t,
    gil::gray32_image_t,
    gil::bgr8_image_t,
    gil::bgr16_image_t,
    gil::bgr32_image_t,
    gil::rgb8_image_t,
    gil::rgb16_image_t,
    gil::rgb32_image_t,
    gil::rgba8_image_t,
    gil::rgba16_image_t,
    gil::rgba32_image_t
>;

template <typename T>
struct consecutive_value
{
    consecutive_value(T start) : current_(start)
    {
        BOOST_TEST(static_cast<int>(current_) >= 0);
    }

    T operator()()
    {
        BOOST_ASSERT(static_cast<int>(current_) + 1 > 0);
        return current_++;
    }

    T current_;
};

template <typename T>
struct reverse_consecutive_value
{
    reverse_consecutive_value(T start) : current_(start)
    {
        BOOST_ASSERT(static_cast<int>(current_) > 0);
    }

    T operator()()
    {
        BOOST_ASSERT(static_cast<int>(current_) + 1 >= 0);
        return current_--;
    }

    T current_;
};

template <typename T>
struct random_value
{
    static_assert(std::is_integral<T>::value, "T must be integral type");
    static constexpr auto range_min = std::numeric_limits<T>::min();
    static constexpr auto range_max = std::numeric_limits<T>::max();

    random_value() : rng_(rd_()), uid_(range_min, range_max) {}

    T operator()()
    {
        auto value = uid_(rng_);
        BOOST_ASSERT(range_min <= value && value <= range_max);
        return static_cast<T>(value);
    }

    std::random_device rd_;
    std::mt19937 rng_;
    std::uniform_int_distribution<typename gil::promote_integral<T>::type> uid_;
};

template <typename Image, typename Generator>
auto generate_image(std::ptrdiff_t size_x, std::ptrdiff_t size_y, Generator&& generate) -> Image
{
    using pixel_t = typename Image::value_type;

    Image out(size_x, size_y);
    gil::for_each_pixel(view(out), [&generate](pixel_t& p) {
        gil::static_generate(p, [&generate]() { return generate(); });
    });

    return out;
}

template <typename Image>
auto create_image(std::ptrdiff_t size_x, std::ptrdiff_t size_y, int channel_value) -> Image
{
    using pixel_t = typename Image::value_type;
    using channel_t = typename gil::channel_type<pixel_t>::type;
    static_assert(std::is_integral<channel_t>::value, "channel must be integral type");

    Image out(size_x, size_y);
    gil::for_each_pixel(view(out), [&channel_value](pixel_t& p) {
        gil::static_fill(p, static_cast<channel_t>(channel_value));
    });

    return out;
}

}}}} // namespace boost::gil::test::fixture

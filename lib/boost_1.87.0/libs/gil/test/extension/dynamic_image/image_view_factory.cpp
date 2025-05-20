//
// Copyright 2022 Marco Langer <langer.m86 at gmail dot com>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//
#include <boost/gil/extension/dynamic_image/algorithm.hpp>
#include <boost/gil/extension/dynamic_image/image_view_factory.hpp>

#include <boost/core/lightweight_test.hpp>

#include <array>

#include "core/image/test_fixture.hpp"
#include "extension/dynamic_image/test_fixture.hpp"

namespace gil = boost::gil;
namespace fixture = boost::gil::test::fixture;

template <std::size_t Channels>
struct generator
{
    generator(int const* data) : data_(data) {}

    auto operator()() -> int
    {
        if(++i_ == Channels) {
            i_ = 0;
            return *data_++;
        }
        else
        {
            return *data_;
        }
    }

    int i_= 0;
    int const* data_;
};

struct test_flipped_up_down_view
{
    template <typename Image>
    void operator()(Image const&)
    {
        using image_t = Image;
        using pixel_t = typename image_t::value_type;
        static constexpr std::size_t num_channels = gil::num_channels<pixel_t>::value;

        std::array<int, 9> pixel_data =
        {
            0, 1, 2,
            3, 4, 5,
            6, 7, 8
        };
        std::array<int, 9> expected_pixel_data =
        {
            6, 7, 8,
            3, 4, 5,
            0, 1, 2
        };

        fixture::dynamic_image source_image =
            fixture::generate_image<image_t>(
                3, 3, generator<num_channels>{pixel_data.data()});

        fixture::dynamic_image expected_image =
            fixture::generate_image<image_t>(
                3, 3, generator<num_channels>{expected_pixel_data.data()});

        auto result_view = gil::flipped_up_down_view(gil::const_view(source_image));

        BOOST_TEST(
            gil::equal_pixels(
                result_view,
                gil::const_view(expected_image)));
    }

    static void run()
    {
        boost::mp11::mp_for_each<fixture::image_types>(test_flipped_up_down_view{});
    }
};

template<typename V1, typename V2>
bool equal_pixels_values(V1&& v1, 
    V2&& v2,
    float threshold = 1e-6f)
{
// convert both images to rgba32f and compare with threshold
   return boost::variant2::visit([=](auto const& v1, auto const& v2) -> bool {
        auto it1 = v1.begin();
        auto it2 = v2.begin();
        while(it1 != v1.end() && it2 != v2.end()) {
            using pixel_t = gil::rgba32f_pixel_t;
            static constexpr std::size_t num_channels = gil::num_channels<pixel_t>::value;
            pixel_t p1{};
            gil::color_convert(*it1++, p1);
            pixel_t p2{};
            gil::color_convert(*it2++, p2);
            for(size_t i = 0; i < num_channels; ++i) {
                if(std::abs(p1[i] - p2[i]) > threshold){
                    return false;
                }
            }
        }
        return true;
    }, 
    std::forward<V1>(v1), 
    std::forward<V2>(v2));
}

struct test_color_converted_view
{
    static void run()
    {
        using dynamic_image_t = gil::any_image<gil::gray8_image_t, gil::gray16_image_t>;
        using src_image_t = gil::gray8_image_t;
        using dst_image_t = gil::gray16_image_t;
        using dst_pixel_t = typename dst_image_t::value_type;
        static constexpr std::size_t num_channels = 1;
        auto color_converter = [](auto const& src, auto& dst) { dst = 2 * src; };
        std::array<int, 9> pixel_data =
        {
            0, 1, 2,
            3, 4, 5,
            6, 7, 8
        };
        std::array<int, 9> expected_pixel_data;
        std::transform(std::begin(pixel_data), 
            std::end(pixel_data), 
            std::begin(expected_pixel_data), 
            [](auto v) { return 2 * v; });

        dynamic_image_t source_image =
        fixture::generate_image<src_image_t>(
                3, 3, generator<num_channels>{pixel_data.data()});

        dynamic_image_t expected_image =
            fixture::generate_image<dst_image_t>(
                3, 3, generator<num_channels>{expected_pixel_data.data()});

        auto result_view = gil::color_converted_view<dst_pixel_t>(gil::const_view(source_image), color_converter);
             
        BOOST_TEST(equal_pixels_values(result_view, 
            gil::const_view(expected_image), 
            1e-6f)); 
    }
};

int main()
{
    test_flipped_up_down_view::run();

    test_color_converted_view::run();

    return ::boost::report_errors();
}

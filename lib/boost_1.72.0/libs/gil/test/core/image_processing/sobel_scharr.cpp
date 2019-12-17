#include <boost/gil/image_processing/numeric.hpp>
#include <boost/gil/detail/math.hpp>
#include <boost/core/lightweight_test.hpp>

#include <algorithm>

namespace gil = boost::gil;

void test_dx_sobel_kernel()
{
    const auto kernel = gil::generate_dx_sobel(1);
    BOOST_TEST(std::equal(kernel.begin(), kernel.end(), gil::dx_sobel.begin()));
}

void test_dx_scharr_kernel()
{
    const auto kernel = gil::generate_dx_scharr(1);
    BOOST_TEST(std::equal(kernel.begin(), kernel.end(), gil::dx_scharr.begin()));
}

void test_dy_sobel_kernel()
{
    const auto kernel = gil::generate_dy_sobel(1);
    BOOST_TEST(std::equal(kernel.begin(), kernel.end(), gil::dy_sobel.begin()));
}

void test_dy_scharr_kernel()
{
    const auto kernel = gil::generate_dy_scharr(1);
    BOOST_TEST(std::equal(kernel.begin(), kernel.end(), gil::dy_scharr.begin()));
}

int main()
{
    test_dx_sobel_kernel();
    test_dx_scharr_kernel();
    test_dy_sobel_kernel();
    test_dy_scharr_kernel();
    return boost::report_errors();
}

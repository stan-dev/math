#ifndef TEST_UNIT_MATH_PRIM_FUN_BETA_REFERENCE_HPP
#define TEST_UNIT_MATH_PRIM_FUN_BETA_REFERENCE_HPP

#include <array>
#include <vector>

namespace stan {
namespace test {
// Reference values of beta(a, b) from mpmath at 50 digits (mp.beta, checked
// against exp(loggamma(a) + loggamma(b) - loggamma(a + b))). Columns: a, b,
// beta(a, b). One argument is large, in both argument orders.
inline const std::vector<std::array<double, 3>> beta_large_reference{{
    {0.5, 1000000.0, 0.0017724540724622611},
    {1000000.0, 0.5, 0.0017724540724622611},
    {0.5, 1000000000.0, 5.6049912170985526e-05},
    {1000000000.0, 0.5, 5.6049912170985526e-05},
    {0.5, 2147483648.0, 3.824811210291909e-05},
    {2147483648.0, 0.5, 3.824811210291909e-05},
    {2.0, 1000000.0, 9.99999000001e-13},
    {1000000.0, 2.0, 9.99999000001e-13},
    {2.0, 1000000000.0, 9.99999999e-19},
    {1000000000.0, 2.0, 9.99999999e-19},
    {2.0, 2147483648.0, 2.168404343961267e-19},
    {2147483648.0, 2.0, 2.168404343961267e-19},
    {10.0, 1000000.0, 3.628636708191183e-55},
    {1000000.0, 10.0, 3.628636708191183e-55},
    {10.0, 1000000000.0, 3.6287998367040043e-85},
    {1000000000.0, 10.0, 3.6287998367040043e-85},
    {10.0, 2147483648.0, 1.7396599603107283e-88},
    {2147483648.0, 10.0, 1.7396599603107283e-88},
}};
}  // namespace test
}  // namespace stan

#endif

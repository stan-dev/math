#ifndef STAN_MATH_PRIM_SCAL_FUN_LGAMMA_STIRLING_DIFF_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LGAMMA_STIRLING_DIFF_HPP

#include <cmath>
//#include <boost/math/special_functions/chebyshev_transform.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/lgamma_stirling.hpp>

namespace stan {
namespace math {

// namespace internal {
//     boost::math::chebyshev_transform<double>
//     build_lgamma_stirling_diff_chebyshev() {
//         auto f = [](double t) {
//             double x = 10 * sqrt(2 / (t + 1));
//             double stirling = 0.5 * log(2*stan::math::pi()) + (x -
//             0.5)*log(x) -x; return (lgamma(x) - stirling);
//         };
//         auto transform = boost::math::chebyshev_transform<double>(f, -1, 1,
//         1e-10); std::cout << "Coeffs: " << transform.coefficients().size() <<
//         std::endl;
//         //std::cout << transform.coefficients()[0];
//         for(const double c : transform.coefficients()) {
//             std::cout << c << " ";
//         }
//         return transform;
//     }

//     boost::math::chebyshev_transform<double> lgamma_stirling_diff_chebyshev =
//     build_lgamma_stirling_diff_chebyshev();
// }

namespace internal {
constexpr double lgamma_stirling_diff_big = 1e5;
}

double lgamma_stirling_diff(const double x) {
  if (x < internal::lgamma_stirling_diff_big) {
    return lgamma(x) - lgamma_stirling(x);
  } else {
    return 1.0 / (12.0 * x) + 1.0 / (288 * x * x);
  }
  // double ten_over_x = 10.0 / x;
  // double t = ten_over_x * ten_over_x * 2 - 1;
  // return internal::lgamma_stirling_diff_chebyshev(t);
}

}  // namespace math
}  // namespace stan

#endif

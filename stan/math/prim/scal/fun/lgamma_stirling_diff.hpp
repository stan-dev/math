#ifndef STAN_MATH_PRIM_SCAL_FUN_LGAMMA_STIRLING_DIFF_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LGAMMA_STIRLING_DIFF_HPP

#include <cmath>
//#include <boost/math/special_functions/chebyshev_transform.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/lgamma_stirling.hpp>


namespace stan {
namespace math {

namespace internal {
    constexpr double lgamma_stirling_diff_big = 1e10;
}


// namespace lgamma_stirling_diff_internal { 
//     typedef long double Real;

//     boost::math::chebyshev_transform<Real> 
//     build_lgamma_stirling_diff_chebyshev() {
//         std::cout << "Size: " << sizeof(Real) << std::endl;
//         auto f = [](Real t) {
//             Real x = 10 * std::sqrt(2 / (t + 1));
//             if(x > internal::lgamma_stirling_diff_big) {
//                 return 1.0 / (12.0 * x);
//             } else {
//                 Real stirling = 0.5 
//                     * std::log(2* static_cast<Real>(stan::math::pi())) 
//                     + (x - 0.5)*std::log(x) -x;
//                 return (std::lgamma(x) - stirling);
//             }
//         };
//         auto transform = boost::math::chebyshev_transform<Real>(f, -1, 1, 1e-8);
//         std::cout << "Coeffs: " << transform.coefficients().size() << std::endl;
//         //std::cout << transform.coefficients()[0];
//         size_t coeffs_to_print = 
//             std::min(size_t(20), transform.coefficients().size());
//         for(size_t i = 0; i < coeffs_to_print ; i ++) {
//             std::cout << std::setprecision(17) << transform.coefficients()[i] 
//                 << std::endl;
//         }
//         return transform;
//     }

//     boost::math::chebyshev_transform<Real> lgamma_stirling_diff_chebyshev = 
//         build_lgamma_stirling_diff_chebyshev();
// }


double lgamma_stirling_diff(const double x) {     
    if (x < internal::lgamma_stirling_diff_big) {
        return std::lgamma(x) - lgamma_stirling(x);
        // double ten_over_x = 10.0 / x;
        // double t = ten_over_x * ten_over_x * 2 - 1;
        // return internal::lgamma_stirling_diff_chebyshev(t); 
    } else {
        return 1.0 / (12.0 * x);
    }
}

}  // namespace math
}  // namespace stan

#endif

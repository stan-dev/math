#ifndef STAN_MATH_PRIM_MAT_FUN_COMMON_TYPE_HPP
#define STAN_MATH_PRIM_MAT_FUN_COMMON_TYPE_HPP

#include <stan/math/prim/arr/fun/common_type.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <boost/math/tools/promotion.hpp>

namespace stan {
  namespace math {
    /**
     * Struct to hold static function for type promotion.
     * This specialization is for matrices of differing types.
     *
     * @tparam T1 type of arg1
     * @tparam T2 type of arg2
     * @tparam R number of rows
     * @tparam C number of columns
     */
    template <typename T1, typename T2, int R, int C>
    struct common_type<Eigen::Matrix<T1, R, C>, Eigen::Matrix<T2, R, C> > {
      typedef Eigen::Matrix<typename common_type<T1, T2>::type, R, C> type;
    };

  }
}

#endif

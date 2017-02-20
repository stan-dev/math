#ifndef STAN_MATH_PRIM_MAT_FUN_PROMOTER_HPP
#define STAN_MATH_PRIM_MAT_FUN_PROMOTER_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
  namespace math {
    /**
     * Struct which holds static function for element-wise type promotion.
     * This specialization promotes Eigen matrix elements of different types.
     *
     * @tparam F type of input element
     * @tparam T type of output element
     * @tparam R number of rows
     * @tparam C number of columns
     */
    template <typename F, typename T, int R, int C>
    struct promoter<Eigen::Matrix<F, R, C>, Eigen::Matrix<T, R, C> > {
      inline static void promote(const Eigen::Matrix<F, R, C>& u,
                          Eigen::Matrix<T, R, C>& t) {
        t.resize(u.rows(), u.cols());
        for (int i = 0; i < u.size(); ++i)
          promoter<F, T>::promote(u(i), t(i));
      }
      inline static Eigen::Matrix<T, R, C>
      promote_to(const Eigen::Matrix<F, R, C>& u) {
        Eigen::Matrix<T, R, C> t;
        promoter<Eigen::Matrix<F, R, C>,
                 Eigen::Matrix<T, R, C> >::promote(u, t);
        return t;
      }
    };

    /**
     * Struct which holds static function for element-wise type promotion.
     * This specialization promotes Eigen matrix elements of same type.
     *
     * @tparam T type of matrix element
     * @tparam R number of rows
     * @tparam C number of columns
     */
    template <typename T, int R, int C>
    struct promoter<Eigen::Matrix<T, R, C>, Eigen::Matrix<T, R, C> > {
      inline static void promote(const Eigen::Matrix<T, R, C>& u,
                          Eigen::Matrix<T, R, C>& t) {
        t = u;
      }
      inline static Eigen::Matrix<T, R, C>
      promote_to(const Eigen::Matrix<T, R, C>& u) {
        return u;
      }
    };

  }
}

#endif

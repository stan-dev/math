#ifndef STAN_MATH_PRIM_MAT_FUN_PROMOTER_HPP
#define STAN_MATH_PRIM_MAT_FUN_PROMOTER_HPP

#include <stan/math/prim/arr/fun/promoter.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
  namespace math {
    /**
     * Struct which holds static function for element-wise type promotion.
     * This specialization promotes Eigen matrix elements of different types.
     *
     * @tparam F type of input element
     * @tparam T type of output element
     */
    template <typename F, typename T>
    struct promoter<Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>,
                    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > {
      inline static void promote(const Eigen::Matrix<F,
                                 Eigen::Dynamic, Eigen::Dynamic>& u,
                                 Eigen::Matrix<T,
                                 Eigen::Dynamic, Eigen::Dynamic>& t) {
        t.resize(u.rows(), u.cols());
        for (int i = 0; i < u.size(); ++i)
          promoter<F, T>::promote(u(i), t(i));
      }
      inline static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
      promote_to(const Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>& u) {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> t;
        promoter<Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>,
                 Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
                 >::promote(u, t);
        return t;
      }
    };

    /**
     * Struct which holds static function for element-wise type promotion.
     * This specialization promotes Eigen matrix elements of same type.
     *
     * @tparam T type of matrix element
     */
    template <typename T>
    struct promoter<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>,
                    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > {
      inline static void promote(const Eigen::Matrix<T,
                                 Eigen::Dynamic, Eigen::Dynamic>& u,
                          Eigen::Matrix<T,
                                 Eigen::Dynamic, Eigen::Dynamic>& t) {
        t = u;
      }
      inline static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
      promote_to(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& u) {
        return u;
      }
    };

  }
}

#endif

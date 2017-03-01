#ifndef STAN_MATH_PRIM_MAT_FUN_PROMOTE_ELEMENTS_HPP
#define STAN_MATH_PRIM_MAT_FUN_PROMOTE_ELEMENTS_HPP

#include <stan/math/prim/arr/fun/promote_elements.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
  namespace math {

    /**
     * Struct with static function for elementwise type promotion.
     *
     * <p>This specialization promotes matrix elements of different types
     * which must be compatible with promotion.
     *
     * @tparam T type of promoted elements
     * @tparam S type of input elements, must be assignable to T
     */
    template <typename T, typename S>
    struct promote_elements<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>,
                            Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> > {
      /**
       * Return input matrix of type S as matrix of type T.
       *
       * @param u matrix of type S, assignable to type T
       * @returns matrix of type T
       */
      inline static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
      promote(const Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>& u) {
        const T* udatap = u.data();
        int C = u.cols();
        int R = u.rows();
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> t(R, C);
        const T* tdatap = t.data();
        for (int i=0, ij=0; i < C; i++)
          for (int j=0; j < R; j++, ij++)
            tdatap[ij] = promote_elements<T, S>::promote(udatap[ij]);
        return t;
      }
    };

    /**
     * Struct with static function for elementwise type promotion.
     *
     * <p>This specialization promotes matrix elements of the same type.
     *
     * @tparam T type of elements
     */
    template <typename T>
    struct promote_elements<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>,
                            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > {
      /**
       * Return input matrix.
       *
       * @param u matrix of type T
       * @returns matrix of type T
       */
      inline static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
      promote(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& u) {
        return u;
      }
    };

  }
}

#endif

#ifndef STAN_MATH_REV_MAT_FUN_MEMALLOC_MATRIX_MAP_HPP
#define STAN_MATH_REV_MAT_FUN_MEMALLOC_MATRIX_MAP_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
  namespace math {

    template <typename T, int R, int C>
    inline Eigen::Map<Eigen::Matrix<T, R, C> >
    memalloc_matrix_map(const Eigen::Matrix<T, R, C>& x) {
      T* x_mem = ChainableStack::memalloc_.template alloc_array<T>(x.size());
      for (int i = 0; i < x.size(); ++i)
        x_mem[i] = x(i);
      return Eigen::Map<Eigen::Matrix<T, R, C> >(x_mem, x.rows(), x.cols());
    }

  }
}
#endif

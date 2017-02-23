#ifndef STAN_MATH_PRIM_MAT_FUN_ARRAY_BUILDER_HPP
#define STAN_MATH_PRIM_MAT_FUN_ARRAY_BUILDER_HPP

#include <stan/math/prim/arr/fun/array_builder.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <vector>

namespace stan {
  namespace math {

    /**
     * Structure for building up arrays in an expression (rather than
     * in statements) using an argumentchaining add() method and
     * a getter method array() to return the result.
     *
     * @tparam T1 type of arg1
     */
    template <typename T, int R, int C>
    struct array_builder<std::vector<Eigen::Matrix<T, R, C> > >{
      std::vector<Eigen::Matrix<T, R, C> >  x_;
      array_builder() : x_() { }
      template <typename F>
      array_builder& add(const Eigen::Matrix<F, R, C>& u) {
        Eigen::Matrix<T, R, C> t;
        promoter<Eigen::Matrix<F, R, C>,
                 Eigen::Matrix<T, R, C> >::promote(u, t);
        x_.push_back(t);
        return *this;
      }
      
      std::vector<Eigen::Matrix<T, R, C> > array() {
        return x_;
      }
    };


  }
}

#endif

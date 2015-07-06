#ifndef STAN_MATH_FWD_MAT_FUN_QUAD_FORM_SYM_HPP
#define STAN_MATH_FWD_MAT_FUN_QUAD_FORM_SYM_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/mat/fun/multiply.hpp>
#include <stan/math/fwd/mat/fun/dot_product.hpp>
#include <stan/math/prim/mat/fun/quad_form_sym.hpp>

namespace stan {
  namespace math {

    template<int RA, int CA, int RB, int CB, typename T>
    inline Eigen::Matrix<stan::math::fvar<T>, CB, CB>
    quad_form_sym(const Eigen::Matrix<stan::math::fvar<T>, RA, CA>& A,
                  const Eigen::Matrix<double, RB, CB>& B) {
      using stan::math::multiply;

      stan::math::check_square("quad_form_sym", "A", A);
      stan::math::check_multiplicable("quad_form_sym",
                                                "A", A,
                                                "B", B);
      stan::math::check_symmetric("quad_form_sym", "A", A);
      Eigen::Matrix<stan::math::fvar<T>, CB, CB> ret(multiply(transpose(B), multiply(A, B)));
      return T(0.5) * (ret + transpose(ret));
    }

    template<int RA, int CA, int RB, typename T>
    inline fvar<T>
    quad_form_sym(const Eigen::Matrix<stan::math::fvar<T>, RA, CA>& A,
                  const Eigen::Matrix<double, RB, 1>& B) {
      using stan::math::multiply;
      using stan::math::dot_product;

      stan::math::check_square("quad_form_sym", "A", A);
      stan::math::check_multiplicable("quad_form_sym",
                                                "A", A,
                                                "B", B);
      stan::math::check_symmetric("quad_form_sym", "A", A);
      return dot_product(B, multiply(A, B));
    }
    template<int RA, int CA, int RB, int CB, typename T>
    inline Eigen::Matrix<stan::math::fvar<T>, CB, CB>
    quad_form_sym(const Eigen::Matrix<double, RA, CA>& A,
                  const Eigen::Matrix<stan::math::fvar<T>, RB, CB>& B) {
      using stan::math::multiply;

      stan::math::check_square("quad_form_sym", "A", A);
      stan::math::check_multiplicable("quad_form_sym",
                                                "A", A,
                                                "B", B);
      stan::math::check_symmetric("quad_form_sym", "A", A);
      Eigen::Matrix<stan::math::fvar<T>, CB, CB> ret(multiply(transpose(B), multiply(A, B)));
      return T(0.5) * (ret + transpose(ret));
    }

    template<int RA, int CA, int RB, typename T>
    inline fvar<T>
    quad_form_sym(const Eigen::Matrix<double, RA, CA>& A,
                  const Eigen::Matrix<stan::math::fvar<T>, RB, 1>& B) {
      using stan::math::multiply;
      using stan::math::dot_product;

      stan::math::check_square("quad_form_sym", "A", A);
      stan::math::check_multiplicable("quad_form_sym",
                                                "A", A,
                                                "B", B);
      stan::math::check_symmetric("quad_form_sym", "A", A);
      return dot_product(B, multiply(A, B));
    }
  }
}

#endif


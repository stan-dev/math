#ifndef STAN_MATH_PRIM_SCAL_META_OPERANDSANDPARTIALS_HPP
#define STAN_MATH_PRIM_SCAL_META_OPERANDSANDPARTIALS_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/meta/VectorView.hpp>

namespace stan {
  namespace math {

    /**
     * This class builds partial derivatives with respect to a set of
     * operands. There are two reason for the generality of this
     * class. The first is to handle vector and scalar arguments
     * without needing to write additional code. The second is to use
     * this class for writing probability distributions that handle
     * primitives, reverse mode, and forward mode variables
     * seamlessly.
     * 
     * The default template class handles the case where the arguments
     * are primitive. There are template specializations for reverse
     * mode and forward mode.
     *
     * @tparam T1 First set of operands.
     * @tparam T2 Second set of operands.
     * @tparam T3 Third set of operands.
     * @tparam T4 Fourth set of operands.
     * @tparam T5 Fifth set of operands.
     * @tparam T6 Sixth set of operands.
     * @tparam T_return_type Return type of the expression. This defaults
     *   to a template metaprogram that calculates the scalar promotion of
     *   T1 -- T6.
     */
    template<typename T1 = double, typename T2 = double, typename T3 = double,
             typename T4 = double, typename T5 = double, typename T6 = double,
             typename T_return_type
             = typename stan::return_type<T1, T2, T3, T4, T5, T6>::type>
    struct OperandsAndPartials {
      VectorView<T_return_type, false, true> d_x1;
      VectorView<T_return_type, false, true> d_x2;
      VectorView<T_return_type, false, true> d_x3;
      VectorView<T_return_type, false, true> d_x4;
      VectorView<T_return_type, false, true> d_x5;
      VectorView<T_return_type, false, true> d_x6;

      /**
       * Constructor.
       *
       * @param x1 first set of operands
       * @param x2 second set of operands
       * @param x3 third set of operands
       * @param x4 fourth set of operands
       * @param x5 fifth set of operands
       * @param x6 sixth set of operands
       */
      OperandsAndPartials(const T1& x1 = 0, const T2& x2 = 0,
                          const T3& x3 = 0, const T4& x4 = 0,
                          const T5& x5 = 0, const T6& x6 = 0) { }

      /**
       * Returns a T_return_type with the value specified with
       * the partial derivatves.
       *
       * @param[in] value Value of the variable
       * @returns a variable with the appropriate value
       */
      T_return_type
      value(double value) {
        return value;
      }
    };

  }
}
#endif

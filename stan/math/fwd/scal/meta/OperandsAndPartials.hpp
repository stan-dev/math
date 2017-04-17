#ifndef STAN_MATH_FWD_SCAL_META_OPERANDSANDPARTIALS_HPP
#define STAN_MATH_FWD_SCAL_META_OPERANDSANDPARTIALS_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/length.hpp>
#include <stan/math/prim/scal/meta/OperandsAndPartials.hpp>
#include <stan/math/fwd/core.hpp>

namespace stan {
  namespace math {

    namespace {
      template <typename T_derivative,
                typename T,
                typename T_partials,
                bool is_vec = is_vector<T>::value,
                bool is_const = is_constant_struct<T>::value>
      struct increment_derivative {
        inline T_derivative operator()(const T& x,
                                       const T_partials& d_dx) {
          return 0;
        }
      };

      template <typename T_derivative,
                typename T,
                typename T_partials>
      struct increment_derivative<T_derivative, T, T_partials, false, false> {
        inline T_derivative operator()(const T& x,
                                       const T_partials& d_dx) {
          return d_dx[0] * x.d_;
        }
      };

      template <typename T_derivative,
                typename T,
                typename T_partials>
      struct increment_derivative<T_derivative, T, T_partials, true, false> {
        inline T_derivative operator()(const T& x,
                                       const T_partials& d_dx) {
          T_derivative derivative(0);
          for (size_t n = 0; n < length(x); n++)
            derivative += d_dx[n] * x[n].d_;
          return derivative;
        }
      };

      template <typename T,
                typename T1, typename D1, typename T2, typename D2,
                typename T3, typename D3, typename T4, typename D4,
                typename T5, typename D5, typename T6, typename D6>
      fvar<T> partials_to_fvar(T& logp,
                               const T1& x1, D1& d_x1,
                               const T2& x2, D2& d_x2,
                               const T3& x3, D3& d_x3,
                               const T4& x4, D4& d_x4,
                               const T5& x5, D5& d_x5,
                               const T6& x6, D6& d_x6) {
        T deriv = 0;
        if (!is_constant_struct<T1>::value)
          deriv += increment_derivative<T, T1, D1>()(x1, d_x1);
        if (!is_constant_struct<T2>::value)
          deriv += increment_derivative<T, T2, D2>()(x2, d_x2);
        if (!is_constant_struct<T3>::value)
          deriv += increment_derivative<T, T3, D3>()(x3, d_x3);
        if (!is_constant_struct<T4>::value)
          deriv += increment_derivative<T, T4, D4>()(x4, d_x4);
        if (!is_constant_struct<T5>::value)
          deriv += increment_derivative<T, T5, D5>()(x5, d_x5);
        if (!is_constant_struct<T6>::value)
          deriv += increment_derivative<T, T6, D6>()(x6, d_x6);
        return fvar<T>(logp, deriv);
      }
    }

    /**
     * This class builds partial derivatives with respect to a set of
     * operands. There are two reason for the generality of this
     * class. The first is to handle vector and scalar arguments
     * without needing to write additional code. The second is to use
     * this class for writing probability distributions that handle
     * primitives, reverse mode, and forward mode variables
     * seamlessly.
     *
     * This is the partial template specialization for when the return
     * type is fvar<T>.
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
    template<typename T1, typename T2, typename T3,
             typename T4, typename T5, typename T6,
             typename T_partials_return>
    struct OperandsAndPartials<T1, T2, T3, T4, T5, T6,
                               fvar<T_partials_return> > {
      typedef fvar<T_partials_return> T_return_type;

      const T1& x1_;
      const T2& x2_;
      const T3& x3_;
      const T4& x4_;
      const T5& x5_;
      const T6& x6_;

      size_t n_partials;
      T_partials_return* all_partials;

      VectorView<T_partials_return,
                 is_vector<T1>::value,
                 is_constant_struct<T1>::value> d_x1;
      VectorView<T_partials_return,
                 is_vector<T2>::value,
                 is_constant_struct<T2>::value> d_x2;
      VectorView<T_partials_return,
                 is_vector<T3>::value,
                 is_constant_struct<T3>::value> d_x3;
      VectorView<T_partials_return,
                 is_vector<T4>::value,
                 is_constant_struct<T4>::value> d_x4;
      VectorView<T_partials_return,
                 is_vector<T5>::value,
                 is_constant_struct<T5>::value> d_x5;
      VectorView<T_partials_return,
                 is_vector<T6>::value,
                 is_constant_struct<T6>::value> d_x6;

      OperandsAndPartials(const T1& x1 = 0, const T2& x2 = 0, const T3& x3 = 0,
                          const T4& x4 = 0, const T5& x5 = 0, const T6& x6 = 0)
        : x1_(x1), x2_(x2), x3_(x3), x4_(x4), x5_(x5), x6_(x6),
          n_partials(!is_constant_struct<T1>::value * length(x1) +
                     !is_constant_struct<T2>::value * length(x2) +
                     !is_constant_struct<T3>::value * length(x3) +
                     !is_constant_struct<T4>::value * length(x4) +
                     !is_constant_struct<T5>::value * length(x5) +
                     !is_constant_struct<T6>::value * length(x6)),
          all_partials(new T_partials_return[n_partials]),
          d_x1(all_partials),
          d_x2(all_partials
               + (!is_constant_struct<T1>::value) * length(x1)),
          d_x3(all_partials
               + (!is_constant_struct<T1>::value) * length(x1)
               + (!is_constant_struct<T2>::value) * length(x2)),
          d_x4(all_partials
               + (!is_constant_struct<T1>::value) * length(x1)
               + (!is_constant_struct<T2>::value) * length(x2)
               + (!is_constant_struct<T3>::value) * length(x3)),
          d_x5(all_partials
               + (!is_constant_struct<T1>::value) * length(x1)
               + (!is_constant_struct<T2>::value) * length(x2)
               + (!is_constant_struct<T3>::value) * length(x3)
               + (!is_constant_struct<T4>::value) * length(x4)),
          d_x6(all_partials
               + (!is_constant_struct<T1>::value) * length(x1)
               + (!is_constant_struct<T2>::value) * length(x2)
               + (!is_constant_struct<T3>::value) * length(x3)
               + (!is_constant_struct<T4>::value) * length(x4)
               + (!is_constant_struct<T5>::value) * length(x5)) {
        std::fill(all_partials, all_partials + n_partials, 0);
      }

      T_return_type
      value(T_partials_return value) {
        return partials_to_fvar(value,
                                x1_, d_x1, x2_, d_x2,
                                x3_, d_x3, x4_, d_x4,
                                x5_, d_x4, x6_, d_x5);
      }

      ~OperandsAndPartials() {
        delete[] all_partials;
      }
    };

  }
}
#endif

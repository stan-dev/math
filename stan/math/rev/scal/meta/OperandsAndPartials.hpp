#ifndef STAN_MATH_REV_SCAL_META_OPERANDSANDPARTIALS_HPP
#define STAN_MATH_REV_SCAL_META_OPERANDSANDPARTIALS_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/length.hpp>
#include <stan/math/prim/scal/meta/OperandsAndPartials.hpp>
#include <stan/math/prim/scal/meta/VectorView.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
  namespace math {

    namespace {
      class partials_vari : public vari {
      private:
        const size_t N_;
        vari** operands_;
        double* partials_;
      public:
        partials_vari(double value,
                      size_t N,
                      vari** operands, double* partials)
          : vari(value),
            N_(N),
            operands_(operands),
            partials_(partials) { }
        void chain() {
          for (size_t n = 0; n < N_; ++n)
            operands_[n]->adj_ += adj_ * partials_[n];
        }
      };

      var partials_to_var(double logp, size_t nvaris,
                          vari** all_varis,
                          double* all_partials) {
        return var(new partials_vari(logp, nvaris, all_varis,
                                     all_partials));
      }

      template<typename T,
               bool is_vec = is_vector<T>::value,
               bool is_const = is_constant_struct<T>::value>
      struct set_varis {
        inline size_t set(vari** /*varis*/, const T& /*x*/) {
          return 0U;
        }
      };
      template<typename T>
      struct set_varis<T, true, false> {
        inline size_t set(vari** varis, const T& x) {
          for (size_t n = 0; n < length(x); n++)
            varis[n] = x[n].vi_;
          return length(x);
        }
      };
      template<>
      struct set_varis<var, false, false> {
        inline size_t set(vari** varis, const var& x) {
          varis[0] = x.vi_;
          return (1);
        }
      };
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
     * type is stan::math::var.
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
             typename T4, typename T5, typename T6>
    struct OperandsAndPartials<T1, T2, T3, T4, T5, T6, stan::math::var> {
      size_t nvaris;
      vari** all_varis;
      double* all_partials;

      VectorView<double,
                 is_vector<T1>::value,
                 is_constant_struct<T1>::value> d_x1;
      VectorView<double,
                 is_vector<T2>::value,
                 is_constant_struct<T2>::value> d_x2;
      VectorView<double,
                 is_vector<T3>::value,
                 is_constant_struct<T3>::value> d_x3;
      VectorView<double,
                 is_vector<T4>::value,
                 is_constant_struct<T4>::value> d_x4;
      VectorView<double,
                 is_vector<T5>::value,
                 is_constant_struct<T5>::value> d_x5;
      VectorView<double,
                 is_vector<T6>::value,
                 is_constant_struct<T6>::value> d_x6;

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
      OperandsAndPartials(const T1& x1 = 0, const T2& x2 = 0, const T3& x3 = 0,
                          const T4& x4 = 0, const T5& x5 = 0, const T6& x6 = 0)
        : nvaris(!is_constant_struct<T1>::value * length(x1) +
                 !is_constant_struct<T2>::value * length(x2) +
                 !is_constant_struct<T3>::value * length(x3) +
                 !is_constant_struct<T4>::value * length(x4) +
                 !is_constant_struct<T5>::value * length(x5) +
                 !is_constant_struct<T6>::value * length(x6)),
          // TODO(carpenter): replace with array allocation fun
          all_varis(static_cast<vari**>
                    (vari::operator new
                     (sizeof(vari*) * nvaris))),
          all_partials(static_cast<double*>
                       (vari::operator new
                        (sizeof(double) * nvaris))),
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
        size_t base = 0;
        if (!is_constant_struct<T1>::value)
          base += set_varis<T1>().set(&all_varis[base], x1);
        if (!is_constant_struct<T2>::value)
          base += set_varis<T2>().set(&all_varis[base], x2);
        if (!is_constant_struct<T3>::value)
          base += set_varis<T3>().set(&all_varis[base], x3);
        if (!is_constant_struct<T4>::value)
          base += set_varis<T4>().set(&all_varis[base], x4);
        if (!is_constant_struct<T5>::value)
          base += set_varis<T5>().set(&all_varis[base], x5);
        if (!is_constant_struct<T6>::value)
          set_varis<T6>().set(&all_varis[base], x6);
        std::fill(all_partials, all_partials+nvaris, 0);
      }

      /**
       * Returns a T_return_type with the value specified with
       * the partial derivatves.
       *
       * @param[in] value Value of the variable
       * @returns a variable with the appropriate value and
       *   the adjoints set for reverse mode autodiff
       */
      stan::math::var value(double value) {
        return partials_to_var(value, nvaris, all_varis,
                               all_partials);
      }
    };

  }
}
#endif

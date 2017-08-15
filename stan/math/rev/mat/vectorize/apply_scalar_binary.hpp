#ifndef STAN_MATH_REV_MAT_VECTORIZE_APPLY_BINARY_SCALAR_HPP
#define STAN_MATH_REV_MAT_VECTORIZE_APPLY_BINARY_SCALAR_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_binary.hpp>
#include <stan/math/rev/core/var.hpp>

namespace stan {

  namespace math {

    /**
     * Template specialization for vectorized functions applying to
     * var and var arguments.
     *
     * @tparam F Type of function defining static apply function.
     */
    template <typename F>
    struct apply_scalar_binary<F, stan::math::var, stan::math::var> {
      /**
       * Function return type, which is <code>var</code>.
       */
      typedef stan::math::var return_t;

      /**
       * Apply the function specified by F to the specified argument.
       * This is defined through a direct application of
       * <code>F::fun()</code>, which must be defined for var and var
       * arguments.
       *
       * @param x var argument scalar.
       * @param y var argument scalar.
       * @return Result of applying F to the scalar x and y.
       */
      static inline return_t apply(const stan::math::var& x,
      const stan::math::var& y) {
        return F::fun(x, y);
      }
    };

    /**
     * Template specialization for vectorized functions applying to
     * var and int arguments.
     *
     * @tparam F Type of function defining static apply function.
     */
    template <typename F>
    struct apply_scalar_binary<F, stan::math::var, int> {
      /**
       * Function return type, which is <code>var</code>.
       */
      typedef stan::math::var return_t;

      /**
       * Apply the function specified by F to the specified argument.
       * This is defined through a direct application of
       * <code>F::fun()</code>, which must be defined for var and int
       * arguments.
       *
       * @param x var argument scalar.
       * @param y int argument scalar.
       * @return Result of applying F to the scalar x and y.
       */
      static inline return_t apply(const stan::math::var& x, int y) {
        return F::fun(x, static_cast<double>(y));
      }
    };

    /**
     * Template specialization for vectorized functions applying to
     * int and var arguments.
     *
     * @tparam F Type of function defining static apply function.
     */
    template <typename F>
    struct apply_scalar_binary<F, int, stan::math::var> {
      /**
       * Function return type, which is <code>var</code>.
       */
      typedef stan::math::var return_t;

      /**
       * Apply the function specified by F to the specified argument.
       * This is defined through a direct application of
       * <code>F::fun()</code>, which must be defined for int and var
       * arguments.
       *
       * @param x int argument scalar.
       * @param y var argument scalar.
       * @return Result of applying F to the scalar x and y.
       */
      static inline return_t apply(int x, const stan::math::var& y) {
        return F::fun(static_cast<double>(x), y);
      }
    };

    /**
     * Template specialization for vectorized functions applying to
     * var and double arguments.
     *
     * @tparam F Type of function defining static apply function.
     */
    template <typename F>
    struct apply_scalar_binary<F, stan::math::var, double> {
      /**
       * Function return type, which is <code>var</code>.
       */
      typedef stan::math::var return_t;

      /**
       * Apply the function specified by F to the specified argument.
       * This is defined through a direct application of
       * <code>F::fun()</code>, which must be defined for var and double
       * arguments.
       *
       * @param x var argument scalar.
       * @param y double argument scalar.
       * @return Result of applying F to the scalar x and y.
       */
      static inline return_t apply(const stan::math::var& x, double y) {
        return F::fun(x, y);
      }
    };

    /**
     * Template specialization for vectorized functions applying to
     * var and double arguments.
     *
     * @tparam F Type of function defining static apply function.
     */
    template <typename F>
    struct apply_scalar_binary<F, double, stan::math::var> {
      /**
       * Function return type, which is <code>var</code>.
       */
      typedef stan::math::var return_t;

      /**
       * Apply the function specified by F to the specified argument.
       * This is defined through a direct application of
       * <code>F::fun()</code>, which must be defined for double and var
       * arguments.
       *
       * @param x double argument scalar.
       * @param y var argument scalar.
       * @return Result of applying F to the scalar x and y.
       */
      static inline return_t apply(double x, const stan::math::var& y) {
        return F::fun(x, y);
      }
    };
  }
}
#endif

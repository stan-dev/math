#ifndef STAN_MATH_REV_MAT_VECTORIZE_APPLY_BINARY_SCALAR_HPP
#define STAN_MATH_REV_MAT_VECTORIZE_APPLY_BINARY_SCALAR_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_binary.hpp>
#include <stan/math/rev/core/var.hpp>

namespace stan {

  namespace math {

    /**
     * Template specialization to var for vectorizing a unary scalar
     * function.  This is a base scalar specialization.  It applies
     * the function specified by the template parameter to the
     * argument.  
     *
     * @tparam F Type of function to apply.
     */
    template <typename F>
    struct apply_scalar_binary<F, stan::math::var, stan::math::var> {
      /**
       * Function return type, which is <code>var</code>.
       */
      typedef stan::math::var return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument variable.
       * @return Function applied to the variable.
       */
      static inline return_t apply(const stan::math::var& x, 
      const stan::math::var& y) {
        return F::fun(x, y);
      }
    };

    /**
     * Template specialization to var for vectorizing a unary scalar
     * function.  This is a base scalar specialization.  It applies
     * the function specified by the template parameter to the
     * argument.  
     *
     * @tparam F Type of function to apply.
     */
    template <typename F>
    struct apply_scalar_binary<F, stan::math::var, int> {
      /**
       * Function return type, which is <code>var</code>.
       */
      typedef stan::math::var return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument variable.
       * @return Function applied to the variable.
       */
      static inline return_t apply(const stan::math::var& x, int y) {
        return F::fun(x, static_cast<double>(y));
      }
    };

    /**
     * Template specialization to var for vectorizing a unary scalar
     * function.  This is a base scalar specialization.  It applies
     * the function specified by the template parameter to the
     * argument.  
     *
     * @tparam F Type of function to apply.
     */
    template <typename F>
    struct apply_scalar_binary<F, int, stan::math::var> {
      /**
       * Function return type, which is <code>var</code>.
       */
      typedef stan::math::var return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument variable.
       * @return Function applied to the variable.
       */
      static inline return_t apply(int x, const stan::math::var& y) {
        return F::fun(static_cast<double>(x), y);
      }
    };

    /**
     * Template specialization to var for vectorizing a unary scalar
     * function.  This is a base scalar specialization.  It applies
     * the function specified by the template parameter to the
     * argument.  
     *
     * @tparam F Type of function to apply.
     */
    template <typename F>
    struct apply_scalar_binary<F, stan::math::var, double> {
      /**
       * Function return type, which is <code>var</code>.
       */
      typedef stan::math::var return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument variable.
       * @return Function applied to the variable.
       */
      static inline return_t apply(const stan::math::var& x, double y) {
        return F::fun(x, y);
      }
    };

    /**
     * Template specialization to var for vectorizing a unary scalar
     * function.  This is a base scalar specialization.  It applies
     * the function specified by the template parameter to the
     * argument.  
     *
     * @tparam F Type of function to apply.
     */
    template <typename F>
    struct apply_scalar_binary<F, double, stan::math::var> {
      /**
       * Function return type, which is <code>var</code>.
       */
      typedef stan::math::var return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument variable.
       * @return Function applied to the variable.
       */
      static inline return_t apply(double x, const stan::math::var& y) {
        return F::fun(x, y);
      }
    };
  }
}
#endif

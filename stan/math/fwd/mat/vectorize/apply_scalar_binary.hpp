#ifndef STAN_MATH_FWD_MAT_VECTORIZE_APPLY_SCALAR_BINARY_HPP
#define STAN_MATH_FWD_MAT_VECTORIZE_APPLY_SCALAR_BINARY_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_binary.hpp>
#include <stan/math/fwd/core/fvar.hpp>

namespace stan {

  namespace math {

    /**
     * Template specialization to fvar<FV> for vectorizing a unary scalar
     * function.  This is a base scalar specialization.  It applies
     * the function specified by the template parameter to the
     * argument.  
     *
     * @tparam F Type of function to apply.
     */
    template <typename F, typename FV>
    struct apply_scalar_binary<F, stan::math::fvar<FV>, 
    stan::math::fvar<FV> > {
      /**
       * Function return type, which is <code>fvar<FV></code>.
       */
      typedef stan::math::fvar<FV> return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument fvar<FV>iable.
       * @return Function applied to the fvar<FV>iable.
       */
      static inline return_t apply(const stan::math::fvar<FV>& x, 
      const stan::math::fvar<FV>& y) {
        return F::fun(x, y);
      }
    };

    /**
     * Template specialization to fvar<FV> for vectorizing a unary scalar
     * function.  This is a base scalar specialization.  It applies
     * the function specified by the template parameter to the
     * argument.  
     *
     * @tparam F Type of function to apply.
     */
    template <typename F, typename FV>
    struct apply_scalar_binary<F, stan::math::fvar<FV>, int> {
      /**
       * Function return type, which is <code>fvar<FV></code>.
       */
      typedef stan::math::fvar<FV> return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument fvar<FV>iable.
       * @return Function applied to the fvar<FV>iable.
       */
      static inline return_t apply(const stan::math::fvar<FV>& x, int y) {
        return F::fun(x, static_cast<double>(y));
      }
    };

    /**
     * Template specialization to fvar<FV> for vectorizing a unary scalar
     * function.  This is a base scalar specialization.  It applies
     * the function specified by the template parameter to the
     * argument.  
     *
     * @tparam F Type of function to apply.
     */
    template <typename F, typename FV>
    struct apply_scalar_binary<F, int, stan::math::fvar<FV> > {
      /**
       * Function return type, which is <code>fvar<FV></code>.
       */
      typedef stan::math::fvar<FV> return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument fvar<FV>iable.
       * @return Function applied to the fvar<FV>iable.
       */
      static inline return_t apply(int x, const stan::math::fvar<FV>& y) {
        return F::fun(static_cast<double>(x), y);
      }
    };

    /**
     * Template specialization to fvar<FV> for vectorizing a unary scalar
     * function.  This is a base scalar specialization.  It applies
     * the function specified by the template parameter to the
     * argument.  
     *
     * @tparam F Type of function to apply.
     */
    template <typename F, typename FV>
    struct apply_scalar_binary<F, stan::math::fvar<FV>, double> {
      /**
       * Function return type, which is <code>fvar<FV></code>.
       */
      typedef stan::math::fvar<FV> return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument fvar<FV>iable.
       * @return Function applied to the fvar<FV>iable.
       */
      static inline return_t apply(const stan::math::fvar<FV>& x, 
      double y) {
        return F::fun(x, y);
      }
    };

    /**
     * Template specialization to fvar<FV> for vectorizing a unary scalar
     * function.  This is a base scalar specialization.  It applies
     * the function specified by the template parameter to the
     * argument.  
     *
     * @tparam F Type of function to apply.
     */
    template <typename F, typename FV>
    struct apply_scalar_binary<F, double, stan::math::fvar<FV> > {
      /**
       * Function return type, which is <code>fvar<FV></code>.
       */
      typedef stan::math::fvar<FV> return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument fvar<FV>iable.
       * @return Function applied to the fvar<FV>iable.
       */
      static inline return_t apply(double x, const 
      stan::math::fvar<FV>& y) {
        return F::fun(x, y);
      }
    };
  }
}
#endif

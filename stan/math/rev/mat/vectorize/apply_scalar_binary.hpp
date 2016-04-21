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
    struct apply_scalar_unary<F, stan::math::var, double> {
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
    struct apply_scalar_unary<F, double, stan::math::var> {
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

    /**
     * Template specialization to var for vectorizing a unary scalar
     * function.  This is a base scalar specialization.  It applies
     * the function specified by the template parameter to the
     * argument.  
     *
     * @tparam F Type of function to apply.
     */
    template <typename F>
    struct apply_scalar_unary<F, stan::math::var, int> {
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
    struct apply_scalar_unary<F, int, stan::math::var> {
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
    template <typename F, typename T>
    struct apply_scalar_unary<F, stan::math::var, T> {
      
      typedef typename Eigen::internal::traits<T>::Scalar scalar_t;

      /**
       * Function return type, which is <code>var</code>.
       */
      typedef Eigen::Matrix<typename apply_scalar_binary<
      F, stan::math::var, scalar_t>::return_t, 
      T::RowsAtCompileTime, T::ColsAtCompileTime> return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument variable.
       * @return Function applied to the variable.
       */
      static inline return_t apply(const stan::math::var& x, const T& y) {
        using stan::math::var;
        return_t result(y.rows(), y.cols());
        for (int j = 0; j < y.cols(); ++j)
          for (int i = 0; i < y.rows(); ++i)
            result(i, j) = apply_scalar_binary<F, var, scalar_t>::
                           apply(x, y(i, j));
        return result;
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
    template <typename F, typename T>
    struct apply_scalar_unary<F, T, stan::math::var> {
      
      typedef typename Eigen::internal::traits<T>::Scalar scalar_t;

      /**
       * Function return type, which is <code>var</code>.
       */
      typedef Eigen::Matrix<typename apply_scalar_binary<
      F, scalar_t, stan::math::var>::return_t, 
      T::RowsAtCompileTime, T::ColsAtCompileTime> return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument variable.
       * @return Function applied to the variable.
       */
      static inline return_t apply(const T& x, const stan::math::var& y) {
        using stan::math::var;
        return_t result(x.rows(), x.cols());
        for (int j = 0; j < x.cols(); ++j)
          for (int i = 0; i < x.rows(); ++i)
            result(i, j) = apply_scalar_binary<F, var, scalar_t>::
                           apply(x(i, j), y);
        return result;
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
    template <typename F, typename T>
    struct apply_scalar_unary<F, stan::math::var, std::vector<T> > {
      /**
       * Function return type, which is <code>var</code>.
       */
      typedef typename std::vector<typename apply_scalar_binary<
      F, stan::math::var, T>::return_t> return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument variable.
       * @return Function applied to the variable.
       */
      static inline return_t apply(const stan::math::var& x, 
      const std::vector<T>& y) {
        using stan::math::var;
        return_t fx(y.size());
        for (size_t i = 0; i < y.size(); ++i)
          fx[i] = apply_scalar_binary<F, var, T>::apply(x, y[i]);
        return fx;
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
    template <typename F, typename T>
    struct apply_scalar_unary<F, std::vector<T>, stan::math::var> {
      /**
       * Function return type, which is <code>var</code>.
       */
      typedef typename std::vector<typename apply_scalar_binary<
      F, T, stan::math::var>::return_t> return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument variable.
       * @return Function applied to the variable.
       */
      static inline return_t apply(const std::vector<T>& x, 
      const stan::math::var& y) {
        using stan::math::var;
        return_t fx(x.size());
        for (size_t i = 0; i < x.size(); ++i)
          fx[i] = apply_scalar_binary<F, var, T>::apply(x[i], y);
        return fx;
      }
    };
  }
}
#endif

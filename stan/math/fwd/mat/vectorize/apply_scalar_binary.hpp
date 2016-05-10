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

    /**
     * Template specialization to fvar<FV> for vectorizing a unary scalar
     * function.  This is a base scalar specialization.  It applies
     * the function specified by the template parameter to the
     * argument.  
     *
     * @tparam F Type of function to apply.
     */
    template <typename F, typename FV, typename T>
    struct apply_scalar_binary<F, stan::math::fvar<FV>, T> {
      
      typedef typename Eigen::internal::traits<T>::Scalar scalar_t;

      /**
       * Function return type, which is <code>fvar<FV></code>.
       */
      typedef Eigen::Matrix<typename apply_scalar_binary<
      F, stan::math::fvar<FV>, scalar_t>::return_t, 
      T::RowsAtCompileTime, T::ColsAtCompileTime> return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument fvar<FV>iable.
       * @return Function applied to the fvar<FV>iable.
       */
      static inline return_t apply(const stan::math::fvar<FV>& x, 
      const T& y) {
        using stan::math::fvar;
        return_t result(y.rows(), y.cols());
        for (int j = 0; j < y.cols(); ++j)
          for (int i = 0; i < y.rows(); ++i)
            result(i, j) = apply_scalar_binary<F, fvar<FV>, scalar_t>::
                           apply(x, y(i, j));
        return result;
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
    template <typename F, typename FV, typename T>
    struct apply_scalar_binary<F, T, stan::math::fvar<FV> > {
      
      typedef typename Eigen::internal::traits<T>::Scalar scalar_t;

      /**
       * Function return type, which is <code>fvar<FV></code>.
       */
      typedef Eigen::Matrix<typename apply_scalar_binary<
      F, scalar_t, stan::math::fvar<FV> >::return_t, 
      T::RowsAtCompileTime, T::ColsAtCompileTime> return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument fvar<FV>iable.
       * @return Function applied to the fvar<FV>iable.
       */
      static inline return_t apply(const T& x, const 
      stan::math::fvar<FV>& y) {
        using stan::math::fvar;
        return_t result(x.rows(), x.cols());
        for (int j = 0; j < x.cols(); ++j)
          for (int i = 0; i < x.rows(); ++i)
            result(i, j) = apply_scalar_binary<F, scalar_t, fvar<FV> >::
                           apply(x(i, j), y);
        return result;
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
    template <typename F, typename FV, typename T>
    struct apply_scalar_binary<F, stan::math::fvar<FV>, std::vector<T> > {
      /**
       * Function return type, which is <code>fvar<FV></code>.
       */
      typedef typename std::vector<typename apply_scalar_binary<
      F, stan::math::fvar<FV>, T>::return_t> return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument fvar<FV>iable.
       * @return Function applied to the fvar<FV>iable.
       */
      static inline return_t apply(const stan::math::fvar<FV>& x, 
      const std::vector<T>& y) {
        using stan::math::fvar;
        return_t fx(y.size());
        for (size_t i = 0; i < y.size(); ++i) {
          fx[i] = apply_scalar_binary<F, fvar<FV>, T>::apply(x, y[i]);
        }
        return fx;
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
    template <typename F, typename T, typename FV>
    struct apply_scalar_binary<F, std::vector<T>, stan::math::fvar<FV> > {
      /**
       * Function return type, which is <code>fvar<FV></code>.
       */
      typedef typename std::vector<typename apply_scalar_binary<
      F, T, stan::math::fvar<FV> >::return_t> return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument fvar<FV>iable.
       * @return Function applied to the fvar<FV>iable.
       */
      static inline return_t apply(const std::vector<T>& x, 
      const stan::math::fvar<FV>& y) {
        using stan::math::fvar;
        return_t fx(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
          fx[i] = apply_scalar_binary<F, T, fvar<FV> >::apply(x[i], y);
        }
        return fx;
      }
    };
  }
}
#endif

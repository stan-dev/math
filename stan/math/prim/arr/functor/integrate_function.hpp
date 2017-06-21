#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_INTEGRATE_1D_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_INTEGRATE_1D_HPP

#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/scal/fun/to_var.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>
#include <stan/math/rev/arr/fun/to_var.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>

#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <boost/bind.hpp>
#include <cmath>
#include <ostream>
#include <vector>
#include <stan/math/prim/arr/functor/DEIntegrator.hpp>

namespace stan {

  namespace math {
    
    template <typename T>
    struct return_type_of_value_of {
      typedef double type;
    };
    
    template <typename T>
    struct return_type_of_value_of <std::vector<T> > {
      typedef std::vector<double> type;
    };
    
    template <typename T, int R, int C>
    struct return_type_of_value_of <Eigen::Matrix<T, R, C> > {
      typedef Eigen::Matrix<double, R, C> type;
    };
    
    template <typename T>
    struct return_type_of_to_var {
      typedef var type;
    };
    
    template <typename T>
    struct return_type_of_to_var <std::vector<T> > {
      typedef std::vector<var> type;
    };
    
    template <typename T, int R, int C>
    struct return_type_of_to_var <Eigen::Matrix<T, R, C> > {
      typedef Eigen::Matrix<var, R, C> type;
    };

    /**
     * Wrap around function to call static method Integrate from
     * class DEIntegrator.
     * 
     * @tparam T Type of f.
     * @param f a functor with signature double (double).
     * @param a lower limit of integration, must be double type.
     * @param b upper limit of integration, must be double type.
     * @param tae target absolute error.
     * @return numeric integral of function f.
     */
    template <typename F>
    inline
    double call_DEIntegrator(const F& f,
                             const double a,
                             const double b,
                             const double tae) {
      return DEIntegrator<F>::Integrate(f, a, b, tae);
    }

    /**
     * Return the numeric integral of a function f given its gradient g.
     * 
     * @tparam T Type of f.
     * @tparam G Type of g.
     * @param f a functor with signature
     * double (double, std::vector<T_beta>) or with signature
     * double (double, T_beta) where the first argument is one being
     * integrated and the second one is either an extra scalar or vector
     * being passed to f.
     * @param g a functor with signature
     * double (double, std::vector<T_beta>, int, std::ostream*) or with
     * signature double (double, T_beta, int, std::ostream*) where the
     * first argument is onebeing integrated and the second one is
     * either an extra scalar or vector being passed to f and the
     * third one selects which component of the gradient vector
     * is to be returned.
     * @param a lower limit of integration, must be double type.
     * @param b upper limit of integration, must be double type.
     * @param msgs stream.
     * @return numeric integral of function f.
     */
    template <typename F, typename G, typename T_beta>
    inline
    typename scalar_type<T_beta>::type
    integrate_1d_grad(const F& f,
                      const G& g,
                      const double a,
                      const double b,
                      const T_beta& beta,
                      std::ostream* msgs) {

      check_finite("integrate_1d", "lower limit", a);
      check_finite("integrate_1d", "upper limit", b);

      //hard case, we want a normalizing factor
      if (!is_constant_struct<T_beta>::value) {
        size_t N = length(beta);
        std::vector<double> grad(N);

        typename return_type_of_value_of<T_beta>::type
          value_of_beta = value_of(beta);
         
        for (size_t n = 0; n < N; n++)
          grad[n] =
            call_DEIntegrator(
              boost::bind<double>(g, _1, value_of_beta,
                                  static_cast<int>(n+1), msgs),
                                  a, b, 1e-6);

        double val_ = call_DEIntegrator(boost::bind<double>(f, _1, value_of_beta,
                                                     msgs),
                                 a, b, 1e-6);

        operands_and_partials<T_beta> ops_partials(beta);
        for (size_t n = 0; n < N; n++)
          ops_partials.edge1_.partials_[n] += grad[n];

        return ops_partials.build(val_);
      //easy case, here we are calculating a normalizing constant,
      //not a normalizing factor, so g doesn't matter at all
      } else {
        return call_DEIntegrator(
          boost::bind<double>(f, _1, value_of(beta), msgs), a, b, 1e-6);
      }
    }

    /**
     * Calculate gradient of f(x, param, std::ostream*)
     * with respect to param_n (which must be an element of param)
     */ 
    template <typename F, typename T_param>
    inline
    double gradient_of_f(const F& f,
                         const double x,
                         const T_param& param,
                         const var& param_n,
                         std::ostream* msgs) {
      set_zero_all_adjoints_nested();
      f(x, param, msgs).grad();
      return param_n.adj();
    }

    /**
     * Return the numeric integral of a function f with its
     * gradients being infered automatically (but slowly).
     * 
     * @tparam T Type of f.
     * @param f a functor with signature
     * double (double, std::vector<T_beta>, std::ostream*) or with
     * signature double (double, T_beta, std::ostream*) where the first
     * argument is one being integrated and the second one is either
     * an extra scalar or vector being passed to f.
     * @param a lower limit of integration, must be double type.
     * @param b upper limit of integration, must be double type.
     * @param msgs stream.
     * @return numeric integral of function f.
     */
    template <typename F, typename T_param>
    inline
    typename scalar_type<T_param>::type integrate_1d(const F& f,
                               const double a,
                               const double b,
                               const T_param& param,
                               std::ostream* msgs) {

      stan::math::check_finite("integrate_1d", "lower limit", a);
      stan::math::check_finite("integrate_1d", "upper limit", b);

      double val_ =
        call_DEIntegrator(
          boost::bind<double>(f,
                              _1, value_of(param), msgs),
                              a, b, 1e-6);

      if (!is_constant_struct<T_param>::value) {
        size_t N = stan::length(param);
        std::vector<double> results(N);

        try {
          start_nested();
          typedef typename return_type_of_to_var<T_param>::type
            clean_T_param;
          clean_T_param clean_param = to_var(value_of(param));
          
          scalar_seq_view<const clean_T_param> clean_param_vec(clean_param);

          for (size_t n = 0; n < N; n++)
            results[n] =
              call_DEIntegrator(boost::bind<double>(gradient_of_f<F, clean_T_param>, f, _1, clean_param, clean_param_vec[n], msgs), a, b, 1e-6);
            
        } catch (const std::exception& e) {
          recover_memory_nested();
          throw;
        }
        recover_memory_nested();

        operands_and_partials<T_param> ops_partials(param);
        for (size_t n = 0; n < N; n++)
          ops_partials.edge1_.partials_[n] += results[n];

        return ops_partials.build(val_);
      } else {
        return val_;
      }
    }

  }

}

#endif

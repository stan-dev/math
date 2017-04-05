#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_INTEGRATE_1D_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_INTEGRATE_1D_HPP

#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/scal/fun/to_var.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>
#include <stan/math/rev/arr/fun/to_var.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>

#include <stan/math/prim/scal/meta/OperandsAndPartials.hpp>
#include <boost/bind.hpp>
#include <cmath>
#include <ostream>
#include <vector>
#include <stan/math/prim/arr/functor/DEIntegrator.hpp>

namespace stan {

  namespace math {

    template <typename T>
    struct return_type_of_value_of {
      typedef T type;
    };
    
    template <>
    struct return_type_of_value_of <var> {
      typedef double type;
    };
    
    template <>
    struct return_type_of_value_of <std::vector<var> > {
      typedef std::vector<double> type;
    };
    
    template <int R, int C>
    struct return_type_of_value_of <Eigen::Matrix<var, R, C> > {
      typedef Eigen::Matrix<double, R, C> type;
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
     * double (double, std::vector<T_beta>, int) or with signature
     * double (double, T_beta, int) where the first argument is one
     * being integrated and the second one is either an extra scalar or
     * vector being passed to f and third one selects which component of
     * the gradient vector is to be returned.
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

        OperandsAndPartials<T_beta> operands_and_partials(beta);
        for (size_t n = 0; n < N; n++)
          operands_and_partials.d_x1[n] += grad[n];

        return operands_and_partials.value(val_);
      //easy case, here we are calculating a normalizing constant,
      //not a normalizing factor, so g doesn't matter at all
      } else {
        return call_DEIntegrator(
          boost::bind<double>(f,
                              _1, value_of(beta), msgs),
                              a, b, 1e-6);
      }
    }


    // ------------------------- integrate_1d ------------------
    // Only functions related to integrate_1d below this point
    
    template <class F, class T_to_var_value_of_beta>
    inline
    double
    partial_f_value_of_beta_n_helper(F& f_, double x, size_t n_,
                                     const T_to_var_value_of_beta
                                       to_var_value_of_beta,
                                     std::ostream* msgs_) {
      f_(x, to_var_value_of_beta, msgs_).grad();
      VectorView<const T_to_var_value_of_beta>
        to_var_value_of_beta_vector_view(to_var_value_of_beta);
      return to_var_value_of_beta_vector_view[n_].adj();
    }

    template <class F, class T_value_of_beta>
    struct partial_f_beta_n {
      const F& f_;
      const T_value_of_beta& value_of_beta_;
      size_t n_;
      std::ostream* msgs_;

      partial_f_beta_n(const F& f,
                       const T_value_of_beta& value_of_beta,
                       std::ostream* msgs)
        : f_(f), value_of_beta_(value_of_beta), msgs_(msgs) {
      }

      void set_n(size_t n) { n_ = n; }

      double operator()(const double x) const {
        start_nested();
        double d_value_of_betan_f;
        try {
          //value_of_beta_ can be a std::vector, an Eigen Vector
          //or a single scalar, so, we call a helper function
          //to deduce the type of to_var(value_of_beta_)
          d_value_of_betan_f = partial_f_value_of_beta_n_helper(f_, x, n_,
                                               to_var(value_of_beta_),
                                               msgs_);
        } catch (const std::exception& /*e*/) {
          recover_memory_nested();
          throw;
        }
        recover_memory_nested();
        return d_value_of_betan_f;
      }

    };


    template <typename F, typename T_value_of_beta>
    inline
    void
    integrate_1d_helper(const F& f,
                              const double a,
                              const double b,
                              const T_value_of_beta&
                                value_of_beta,
                              size_t N,
                              std::vector<double>& grad,
                              std::ostream* msgs) {

        partial_f_beta_n<F, const T_value_of_beta>
          instance_of_partial_f_beta_n(f, value_of_beta, msgs);

        for (size_t n = 0; n < N; n++) {
          instance_of_partial_f_beta_n.set_n(n);
          grad[n] = call_DEIntegrator(instance_of_partial_f_beta_n,
                                          a, b, 1e-6);
        }

    }


    /**
     * Return the numeric integral of a function f with its
     * gradients being infered automatically (but slowly).
     * 
     * @tparam T Type of f.
     * @param f a functor with signature
     * double (double, std::vector<T_beta>) or with signature
     * double (double, T_beta) where the first argument is one being
     * integrated and the second one is either an extra scalar or vector
     * being passed to f.
     * @param a lower limit of integration, must be double type.
     * @param b upper limit of integration, must be double type.
     * @param msgs stream.
     * @return numeric integral of function f.
     */
    template <typename F, typename T_beta>
    inline
    typename scalar_type<T_beta>::type
    integrate_1d(const F& f,
                       const double a,
                       const double b,
                       const T_beta& beta,
                       std::ostream* msgs) {

      check_finite("integrate_1d", "lower limit", a);
      check_finite("integrate_1d", "upper limit", b);


      double val_ =
        call_DEIntegrator(boost::bind<double>(f,
                                                  _1,
                                                  value_of(beta),
                                                  msgs),
                              a, b, 1e-6);

      if (!is_constant_struct<T_beta>::value) {
        size_t N = length(beta);
        std::vector<double> grad(N);

        //beta can be a std::vector, an Eigen Vector
        //or a single scalar, so, we call a helper function
        //to deduce the type of value_of(beta)
        integrate_1d_helper(f, a, b, value_of(beta),
                                  N, grad, msgs);

        OperandsAndPartials<T_beta> operands_and_partials(beta);
        for (size_t n = 0; n < N; n++) {
          operands_and_partials.d_x1[n] += grad[n];
        }

        return operands_and_partials.value(val_);
      } else {
        return val_;
      }
      
    }



  }

}

#endif

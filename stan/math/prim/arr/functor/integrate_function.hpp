#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_INTEGRATE_FUNCTION_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_INTEGRATE_FUNCTION_HPP

#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>

#include <stan/math/prim/scal/meta/OperandsAndPartials.hpp>
#include <boost/bind.hpp>
#include <cmath>
#include <boost/lambda/lambda.hpp>
#include <ostream>
#include <vector>
#include <stan/math/prim/arr/functor/DEIntegrator.hpp>

namespace stan {

  namespace math {

    namespace {

      //to_var functions: construct a clean
      //(vector/matrix of) vars from a (vector/matrix of) doubles
      template <int R, int C>
      inline typename Eigen::Matrix<var, R, C>
      to_var(const Eigen::Matrix<double, R, C>& x) {
        int S = x.size();
        Eigen::Matrix<var, R, C> result(x.rows(), x.cols());
        var* datap = result.data();
        const double* datax = x.data();
        for (int i=0; i < S; i++)
          datap[i] = var(datax[i]);
        return result;
      }

      inline
      std::vector<var> to_var(const std::vector<double>& x) {
        size_t size = x.size();
        std::vector<var> result(size);
        for (int i=0; i < size; i++)
          result[i] = var(x[i]);
        return result;
      }

      inline var to_var(const double x) {
        return var(x);
      }

      template <typename G>
      inline
      double integrate_definite_1d(const G& g,
                                   const double a,
                                   const double b,
                                   const double tae) {
        return DEIntegrator<G>::Integrate(g, a, b, tae);
      }

    }

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
    integrate_function_helper(const F& f,
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
          grad[n] = integrate_definite_1d(instance_of_partial_f_beta_n,
                                          a, b, 1e-6);
        }

    }


    template <typename F, typename T_beta>
    inline
    typename scalar_type<T_beta>::type
    integrate_function(const F& f,
                       const double a,
                       const double b,
                       const T_beta& beta,
                       std::ostream* msgs) {

      check_finite("integrate_function", "lower limit", a);
      check_finite("integrate_function", "upper limit", b);


      double val_ =
        integrate_definite_1d(boost::bind<double>(f,
                                                  boost::lambda::_1,
                                                  value_of(beta),
                                                  msgs),
                              a, b, 1e-6);

      if (!is_constant_struct<T_beta>::value) {
        size_t N = length(beta);
        std::vector<double> grad(N);

        //beta can be a std::vector, an Eigen Vector
        //or a single scalar, so, we call a helper function
        //to deduce the type of value_of(beta)
        integrate_function_helper(f, a, b, value_of(beta),
                                  N, grad, msgs);

        OperandsAndPartials<T_beta> operands_and_partials(beta);
        for (size_t n = 0; n < N; n++) {
          operands_and_partials.d_x1[n] += grad[n];
        }

        return operands_and_partials.value(val_);
      } else
        return val_;
    }



    template <typename F, typename G, typename T_value_of_beta>
    inline
    double integrate_function_grad_helper(const F& f,
                                          const G& g,
                                          const double a,
                                          const double b,
                                          const T_value_of_beta&
                                            value_of_beta,
                                          const size_t N,
                                          std::vector<double>& grad,
                                          std::ostream* msgs) {

      for (size_t n = 0; n < N; n++)
        grad[n] =
        integrate_definite_1d(
          boost::bind<double>(g,
                              boost::lambda::_1, value_of_beta,
                              static_cast<int>(n+1), msgs),
                              a, b, 1e-6);

      return integrate_definite_1d(
               boost::bind<double>(f,
                                   boost::lambda::_1, value_of_beta, msgs),
                                   a, b, 1e-6);
    }

    template <typename F, typename G, typename T_beta>
    inline
    typename scalar_type<T_beta>::type
    integrate_function_grad(const F& f,
                            const G& g,
                            const double a,
                            const double b,
                            const T_beta& beta,
                            std::ostream* msgs) {

      check_finite("integrate_function", "lower limit", a);
      check_finite("integrate_function", "upper limit", b);

      if (!is_constant_struct<T_beta>::value) {
        size_t N = length(beta);
        std::vector<double> grad(N);

        double val_ =
          integrate_function_grad_helper(f, g, a, b, value_of(beta),
                                         N, grad, msgs);


        OperandsAndPartials<T_beta> operands_and_partials(beta);
        for (size_t n = 0; n < N; n++)
          operands_and_partials.d_x1[n] += grad[n];

        return operands_and_partials.value(val_);
      } else
        return integrate_definite_1d(
          boost::bind<double>(f,
                              boost::lambda::_1, value_of(beta), msgs),
                              a, b, 1e-6);
    }


  }

}

#endif

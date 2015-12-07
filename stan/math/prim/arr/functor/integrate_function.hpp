#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_INTEGRATE_FUNCTION_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_INTEGRATE_FUNCTION_HPP

#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_nonzero_size.hpp>
#include <stan/math/prim/mat/err/check_ordered.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>

#include <stan/math/prim/scal/meta/OperandsAndPartials.hpp>
#include <stan/math/prim/scal/meta/contains_nonconstant_struct.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/scal/fun/inc_beta.hpp>
#include <boost/bind.hpp>
#include <cmath>
#include <boost/lambda/lambda.hpp>
#include <ostream>
#include <vector>
#include <stan/math/prim/arr/functor/DEIntegrator.hpp>

namespace stan {

  namespace math {

    namespace {

      template <int R, int C>
      inline typename Eigen::Matrix<var, R, C>
      to_var(const Eigen::Matrix<double, R, C>& x) {
        int R_ = x.rows();
        int C_ = x.cols();
        int S = x.size();
        Eigen::Matrix<var, R, C> result(R_, C_);
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

      template <typename F>
      inline
      double get_integrate_val(const F& f,
                               const double a,
                               const double b,
                               const double tae) {
        return DEIntegrator<F>::Integrate(f, a, b, tae);
      }

      template <typename F, typename T_param>
      inline
      double get_integrate_gradient(const F& f,
                                    const double x,
                                    const T_param& param,
                                    const var& param_n,
                                    std::ostream* msgs) {

        set_zero_all_adjoints_nested();
        f(x, param, msgs).grad();
        return param_n.adj();
      }

      template <typename F, typename T_param>
      inline
      void do_nested_grad_call(const F& f,
                                 const double a,
                                 const double b,
                                 const T_param& param,
                                 const size_t N,
                                 std::vector<double>& results,
                                 std::ostream* msgs) {

        VectorView<const T_param> param_vec(param);

        for (size_t n = 0; n < N; n++)
          results[n] =
          get_integrate_val(boost::bind<double>(get_integrate_gradient<F, T_param>, f, boost::lambda::_1, param, param_vec[n], msgs), a, b, 1e-6);

      }

    }

    template <typename F, typename T_param>
    inline
    typename scalar_type<T_param>::type integrate_function(const F& f,
                               const double a,
                               const double b,
                               const T_param& param,
                               std::ostream* msgs) {

      stan::math::check_finite("integrate_function", "lower limit", a);
      stan::math::check_finite("integrate_function", "upper limit", b);


      //FIXME: convert vector of var to vector of double
      double val_ =
        get_integrate_val(
          boost::bind<double>(f,
                              boost::lambda::_1, value_of(param), msgs),
                              a, b, 1e-6);

      if (!is_constant_struct<T_param>::value) {

        size_t N = stan::length(param);
        std::vector<double> results(N);

        try {
          start_nested();
          do_nested_grad_call(f, a, b, to_var(value_of(param)),
          N, results, msgs);
        } catch (const std::exception& e) {
          recover_memory_nested();
          throw;
        }
        recover_memory_nested();

        OperandsAndPartials<T_param> operands_and_partials(param);
        for (size_t n = 0; n < N; n++)
          operands_and_partials.d_x1[n] += results[n];

        return operands_and_partials.to_var(val_, param);
      } else
        return val_;
    }

  }

}

#endif

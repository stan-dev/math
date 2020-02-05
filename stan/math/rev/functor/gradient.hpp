#ifndef STAN_MATH_REV_FUNCTOR_GRADIENT_HPP
#define STAN_MATH_REV_FUNCTOR_GRADIENT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stdexcept>

namespace stan {
namespace math {

/**
 * Calculate the value and the gradient of the specified function
 * at the specified argument.
 *
 * <p>The functor must implement
 *
 * <code>
 * var
 * operator()(const
 * Eigen::Matrix<var, Eigen::Dynamic, 1>&)
 * </code>
 *
 * using only operations that are defined for
 * <code>var</code>.  This latter constraint usually
 * requires the functions to be defined in terms of the libraries
 * defined in Stan or in terms of functions with appropriately
 * general namespace imports that eventually depend on functions
 * defined in Stan.
 *
 * <p>Time and memory usage is on the order of the size of the
 * fully unfolded expression for the function applied to the
 * argument, independently of dimension.
 *
 * @tparam F Type of function
 * @param[in] f Function
 * @param[in] x Argument to function
 * @param[out] fx Function applied to argument
 * @param[out] grad_fx Gradient of function at argument
 */
template <typename F>
void gradient(const F& f, const Eigen::Matrix<double, Eigen::Dynamic, 1>& x,
              double& fx, Eigen::Matrix<double, Eigen::Dynamic, 1>& grad_fx) {
  start_nested();
  try {
    Eigen::Matrix<var, Eigen::Dynamic, 1> x_var(x.size());
    for (int i = 0; i < x.size(); ++i)
      x_var(i) = var(new vari(x(i), false));
    var fx_var = f(x_var);
    fx = fx_var.val();
    grad_fx.resize(x.size());
    grad(fx_var.vi_);
    grad_fx = x_var.adj();
  } catch (const std::exception& /*e*/) {
    recover_memory_nested();
    throw;
  }
  recover_memory_nested();
}

template<size_t N, typename F, typename... Ts>
class invoke_packed {
  public:
  static var invoke(const F& f, const Eigen::Matrix<var, Eigen::Dynamic, 1>& x_arr, Ts... xs) {
    return invoke_packed<N - 1, F, const var&, Ts...>::invoke(f, x_arr, x_arr(x_arr.size() - N), xs...);
  }
};

template<typename F, typename... Ts>
class invoke_packed<0, F, Ts...> {
  public:
  static var invoke(const F& f, const Eigen::Matrix<var, Eigen::Dynamic, 1>& , Ts... xs) {
    return f(xs...);
  }
};


template <size_t N, typename F> 
class param_packed {
  public:
    param_packed(const F& f) : f_(f) {};

    var operator() (const Eigen::Matrix<var, Eigen::Dynamic, 1>& x) const {
      // auto f_bound = std::bind(f_, x(x.size() - N - 1));
      // param_packed<N - 1, decltype(f_bound)> p(f_bound);
      // return p(x);
      return invoke_packed<N,F>::invoke(f_, x);
    }

  private:
    const F& f_;
};

// template <typename F> 
// class param_packed<0, F> {
//   public:
//     param_packed(const F& f) : f_(f) {};

//     var operator() (const Eigen::Matrix<var, Eigen::Dynamic, 1>& x) const {
//       return f_();
//     }
//   private:
//     const F& f_;
// };


template <typename F, typename... Ts>
void gradient(const F& f, 
              double& fx, Eigen::Matrix<double, Eigen::Dynamic, 1>& grad_fx,
              Ts... xs) {
  //double x_vec[] = {xs...};
  Eigen::Matrix<double, sizeof...(xs), 1> x(xs...);
  // x.resize(sizeof...(xs));
  // x << x_vec;        
  gradient(param_packed<sizeof...(xs), F>(f), x, fx, grad_fx);
}

}  // namespace math
}  // namespace stan
#endif

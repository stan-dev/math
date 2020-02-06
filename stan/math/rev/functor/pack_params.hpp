#ifndef STAN_MATH_REV_FUNCTOR_PACK_PARAMS_HPP
#define STAN_MATH_REV_FUNCTOR_PACK_PARAMS_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stdexcept>

namespace stan {
namespace math {

template <size_t N, typename F, typename... Ts>
class invoke_packed {
 public:
  static var invoke(const F& f,
                    const Eigen::Matrix<var, Eigen::Dynamic, 1>& x_arr,
                    Ts... xs) {
    return invoke_packed<N - 1, F, const var&, Ts...>::invoke(
        f, x_arr, x_arr(x_arr.size() - N), xs...);
  }
};

template <typename F, typename... Ts>
class invoke_packed<0, F, Ts...> {
 public:
  static var invoke(const F& f, const Eigen::Matrix<var, Eigen::Dynamic, 1>&,
                    Ts... xs) {
    return f(xs...);
  }
};

template <size_t N, typename F>
class pack_params {
 public:
  explicit pack_params(const F& f) : f_(f) {}

  var operator()(const Eigen::Matrix<var, Eigen::Dynamic, 1>& x) const {
    return invoke_packed<N, F>::invoke(f_, x);
  }

 private:
  const F& f_;
};

}
}
#endif

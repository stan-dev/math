#ifndef STAN_MATH_REV_FUN_EXP_HPP
#define STAN_MATH_REV_FUN_EXP_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/cos.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/rev/fun/is_inf.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/isfinite.hpp>
#include <stan/math/rev/fun/is_nan.hpp>
#include <stan/math/rev/fun/sin.hpp>
#include <cmath>
#include <complex>
#include <tbb/task_arena.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace stan {
namespace math {

namespace internal {
class exp_vari : public op_v_vari {
 public:
  explicit exp_vari(vari* avi) : op_v_vari(std::exp(avi->val_), avi) {}
  void chain() { avi_->adj_ += adj_ * val_; }
};
}  // namespace internal

/**
 * Return the exponentiation of the specified variable (cmath).
 *
   \f[
   \mbox{exp}(x) =
   \begin{cases}
     e^x & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{exp}(x)}{\partial x} =
   \begin{cases}
     e^x & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable to exponentiate.
 * @return Exponentiated variable.
 */
inline var exp(const var& a) { return var(new internal::exp_vari(a.vi_)); }

template <typename F, typename T>
inline auto app_func(const F& f, T&& x) {
    // Run nested autodiff in this scope
    Eigen::MatrixXd res(2, x.rows());

    tbb::parallel_for(
      tbb::blocked_range<size_t>(0, x.size()), 
      [&x,&res,&f](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i < r.end(); ++i) {
          nested_rev_autodiff nested;
          var in = deep_copy_vars(x[i]);
          var out = f(in);
          out.grad();
          res(0, i) = out.val();
          res(1, i) = in.adj();
        }
      });
  return res;
}

template <typename Container,
          require_eigen_col_vector_st<is_var, Container>* = nullptr>
inline auto exp(Container&& x) {
  Eigen::Matrix<var,-1,1> result(x.rows());
  auto f = [&](const auto& xi) { return stan::math::exp(xi); };
  auto out = app_func(f, std::forward<Container>(x));

  for(int i = 0; i < x.rows(); ++i) {
    result[i] = var(new precomp_v_vari(out(0, i), x[i].vi_, out(1, i)));
  }

  return result;
}

/**
 * Return the exponentiation (base e) of the specified complex number.
 * @param z argument
 * @return exponentiation of argument
 */
inline std::complex<var> exp(const std::complex<var>& z) {
  return internal::complex_exp(z);
}

}  // namespace math
}  // namespace stan
#endif

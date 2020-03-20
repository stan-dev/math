#ifndef STAN_MATH_FWD_FUN_COS_HPP
#define STAN_MATH_FWD_FUN_COS_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the cosine of a variable.
 *
 * The derivative is defined by
 *
   \f[
   \frac{\partial\, \mbox{cos}(x)}{\partial x} =
   \begin{cases}
     -\sin(x) & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T Inner type of the fvar
 * @param x fvar<T> variable
 * @return Cosine of variable.
 */
template <typename T>
inline fvar<T> cos(const fvar<T>& x) {
  using std::cos;
  using std::sin;
  return fvar<T>(cos(x.val_), x.d_ * -sin(x.val_));
}

/**
 * Return the cosine of each variable in a container.
 *
 * @tparam Container Type of container
 * @param x Container of fvar
 * @return Cosine of each variable in container.
 */
template <typename Container,
          require_container_st<is_container, is_fvar, Container>...>
inline auto cos(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) {
        using T_plain = plain_type_t<decltype(v)>;
        const Eigen::Ref<const T_plain>& v_ref = v;
        auto vals = v_ref.val().eval();

        T_plain result(v_ref.rows(), v_ref.cols());
        result.val() = cos(vals);
        result.d().array() = v_ref.d().array() * -sin(vals).array();

        return result;
});
}
}  // namespace math
}  // namespace stan
#endif

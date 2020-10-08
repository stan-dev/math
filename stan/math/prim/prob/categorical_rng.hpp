#ifndef STAN_MATH_PRIM_PROB_CATEGORICAL_RNG_HPP
#define STAN_MATH_PRIM_PROB_CATEGORICAL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/cumulative_sum.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

template <typename T_t, class RNG>
inline typename VectorBuilder<true, int, value_type_t<T_t>>::type categorical_rng(
    const T_t& t, RNG& rng) {
  using boost::uniform_01;
  using boost::variate_generator;

  vector_seq_view<T_t> t_vec(t);
  size_t N = stan::math::size_mvt(t);
  
  Eigen::VectorXd index(t_vec[0].size());
  VectorBuilder<true, int, value_type_t<T_t>> output(N);
  variate_generator<RNG&, uniform_01<> > uniform01_rng(rng, uniform_01<>());

  for (size_t n = 0; n < N; ++n) {
    check_simplex("categorical_rng", "Probabilities parameter", t_vec[n]);

    index.setZero();
    index = cumulative_sum(t_vec[n]);

    double c = uniform01_rng();
    int b = 0;
    while (c > index(b)) {
      b++;
    }
    output[n] = b + 1;
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif

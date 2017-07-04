#ifndef STAN_MATH_PRIM_SCAL_PROB_NORMAL_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NORMAL_RNG_HPP

#include <stan/math/prim/mat/meta/promote_vector.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view_writable.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/max_size.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>

namespace stan {
  namespace math {

    template <typename T_loc, typename T_scale, class RNG>
    inline typename promote_vector<T_loc, T_scale>::type
    normal_rng(const T_loc &mu,
               const T_scale &sigma,
               RNG& rng) {
      using boost::variate_generator;
      using boost::normal_distribution;

      static const char* function("normal_rng");

      scalar_seq_view<T_loc> mu_vec(mu);
      scalar_seq_view<T_scale> sigma_vec(sigma);

      check_finite(function, "Location parameter", mu);
      check_not_nan(function, "Location parameter", mu);
      check_positive(function, "Scale parameter", sigma);
      check_not_nan(function, "Scale parameter", sigma);
      if (mu_vec.size() > 1 && sigma_vec.size() > 1)
        check_size_match(function,
                         "Location parameter", mu_vec.size(),
                         "Scale Parameter", sigma_vec.size());

      size_t N = max_size(mu, sigma);

      typename promote_vector<T_loc, T_scale>::type output(N);
      scalar_seq_view_writable<typename promote_vector<T_loc, T_scale>::type>
        output_writer(output);

      for (size_t n = 0; n < N; n++) {
        variate_generator<RNG&, normal_distribution<> >
          norm_rng(rng, normal_distribution<>(mu_vec[n], sigma_vec[n]));
        output_writer[n] = norm_rng();
      }

      return output;
    }
  }
}
#endif

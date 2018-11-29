#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_COUPLED_ODE_OBSERVER_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_COUPLED_ODE_OBSERVER_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/broadcast_array.hpp>

// grr... need to pull in the MAT operands_and_partials to make things
// work with vector... should be moved!
#include <stan/math/rev/mat/meta/operands_and_partials.hpp>
#include <stan/math/prim/arr/fun/sum.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

#include <stan/math/prim/scal/meta/return_type.hpp>

#include <vector>

namespace stan {
namespace math {

/**
 * Observer for the coupled states.  Holds a reference to
 * an externally defined vector of vectors passed in at
 * construction time.
 */
template <typename T1, typename T2>
struct coupled_ode_observer {
  typedef typename stan::return_type<T1, T2>::type return_t;

  typedef operands_and_partials<std::vector<T1>, std::vector<T2>>
      ops_partials_t;

  std::vector<std::vector<return_t>>& y_;
  // const std::vector<T1>& y0_;
  // const std::vector<T2>& theta_;
  const std::size_t num_outputs_;
  const std::size_t N_;
  const std::size_t M_;
  const std::size_t index_offset_theta_;
  std::size_t n_;
  ops_partials_t ops_partials_;

  /**
   * Construct a coupled ODE observer from the specified coupled
   * vector.
   *
   * @param y_coupled reference to a vector of vector of doubles.
   */
  explicit coupled_ode_observer(std::vector<std::vector<return_t>>& y,
                                const std::vector<T1>& y0,
                                const std::vector<T2>& theta,
                                std::size_t num_outputs)
      : y_(y),
        // y0_(y0), theta_(theta),
        num_outputs_(num_outputs),
        N_(y0.size()),
        M_(theta.size()),
        index_offset_theta_(is_constant_struct<T1>::value ? 0 : N_ * N_),
        n_(0),
        ops_partials_(y0, theta) {}

  /**
   * Callback function for Boost's ODE solver to record values.
   *
   * @param coupled_state solution at the specified time.
   * @param t time of solution.
   */
  void operator()(const std::vector<double>& coupled_state, double t) {
    if (n_++ == 0)
      return;

    std::vector<return_t> yt;
    yt.reserve(N_);

    for (size_t j = 0; j < N_; j++) {
      // iterate over parameters for each equation
      if (!is_constant_struct<T1>::value) {
        for (std::size_t k = 0; k < N_; k++)
          ops_partials_.edge1_.partials_[k] = coupled_state[N_ + N_ * k + j];
      }

      if (!is_constant_struct<T2>::value) {
        for (std::size_t k = 0; k < M_; k++)
          ops_partials_.edge2_.partials_[k]
              = coupled_state[N_ + index_offset_theta_ + N_ * k + j];
      }

      yt.emplace_back(ops_partials_.build(coupled_state[j]));
    }

    y_.emplace_back(yt);
  }
};

}  // namespace math

}  // namespace stan

#endif

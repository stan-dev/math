#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_COUPLED_ODE_OBSERVER_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_COUPLED_ODE_OBSERVER_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/broadcast_array.hpp>

// grr... need to pull in the MAT operands_and_partials to make things
// work with vector... should be moved!
#include <stan/math/rev/mat/meta/operands_and_partials.hpp>

#include <stan/math/rev/scal/meta/operands_and_partials.hpp>
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
template <typename F, typename T1, typename T2, typename T_t0, typename T_ts>
struct coupled_ode_observer {
  typedef typename stan::return_type<T1, T2, T_t0, T_ts>::type return_t;

  typedef operands_and_partials<std::vector<T1>, std::vector<T2>, T_t0, T_ts>
      ops_partials_t;

  const F& f_;
  const std::vector<T1>& y0_;
  const T_t0& t0_;
  const std::vector<T_ts>& ts_;
  const std::vector<T2>& theta_;
  const std::vector<double>& x_;
  const std::vector<int>& x_int_;
  std::ostream* msgs_;
  std::vector<std::vector<return_t>>& y_;
  const std::size_t N_;
  const std::size_t M_;
  const std::size_t index_offset_theta_;
  std::size_t n_;

  /**
   * Construct a coupled ODE observer from the specified coupled
   * vector.
   *
   * @param y_coupled reference to a vector of vector of doubles.
   */
  explicit coupled_ode_observer(const F& f, const std::vector<T1>& y0,
                                const std::vector<T2>& theta, const T_t0& t0,
                                const std::vector<T_ts>& ts,
                                const std::vector<double>& x,
                                const std::vector<int>& x_int,
                                std::ostream* msgs,
                                std::vector<std::vector<return_t>>& y)
      : f_(f),
        y0_(y0),
        t0_(t0),
        ts_(ts),
        theta_(theta),
        x_(x),
        x_int_(x_int),
        msgs_(msgs),
        y_(y),
        N_(y0.size()),
        M_(theta.size()),
        index_offset_theta_(is_constant_struct<T1>::value ? 0 : N_ * N_),
        n_(0) {}

  /**
   * Callback function for ODE solvers to record values.
   *
   * @param coupled_state solution at the specified time.
   * @param t time of solution.
   */
  void operator()(const std::vector<double>& coupled_state, double t) {
    if (n_++ == 0)
      return;

    std::vector<return_t> yt;
    yt.reserve(N_);

    ops_partials_t ops_partials(y0_, theta_, t0_, ts_[n_ - 2]);

    std::vector<double> dy_dt;
    if (!is_constant_struct<T_ts>::value) {
      std::vector<double> y_dbl(coupled_state.begin(),
                                coupled_state.begin() + N_);
      dy_dt = f_(value_of(ts_[n_ - 2]), y_dbl, value_of(theta_), x_, x_int_,
                 msgs_);
      check_size_match("coupled_ode_observer", "dy_dt", dy_dt.size(), "states",
                       N_);
    }

    for (size_t j = 0; j < N_; j++) {
      // iterate over parameters for each equation
      if (!is_constant_struct<T1>::value) {
        for (std::size_t k = 0; k < N_; k++)
          ops_partials.edge1_.partials_[k] = coupled_state[N_ + N_ * k + j];
      }

      if (!is_constant_struct<T2>::value) {
        for (std::size_t k = 0; k < M_; k++)
          ops_partials.edge2_.partials_[k]
              = coupled_state[N_ + index_offset_theta_ + N_ * k + j];
      }

      if (!is_constant_struct<T_ts>::value) {
        ops_partials.edge4_.partials_[0] = dy_dt[j];
      }

      yt.emplace_back(ops_partials.build(coupled_state[j]));
    }

    y_.emplace_back(yt);
  }
};

}  // namespace math

}  // namespace stan

#endif

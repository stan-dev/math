#ifndef STAN_MATH_PRIM_SCAL_FUN_CHOOSE_HPP
#define STAN_MATH_PRIM_SCAL_FUN_CHOOSE_HPP

#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <limits>

namespace stan {
  namespace math {
    /**
     * Return the binomial coefficient for the specified integer
     * arguments.
     *
     * The binomial coefficient, \f${n \choose k}\f$, read "n choose k", is
     * defined for \f$0 \leq k \leq n\f$ (otherwise return 0) by
     *
     * \f${n \choose k} = \frac{n!}{k! (n-k)!}\f$.
     *
     * @param n total number of objects.
     * @param k number of objects chosen.
     * @return n choose k or 0 iff k > n unless it overflows
     */
    inline int
    choose(const int n, const int k) {
      using stan::math::check_nonnegative;
      using stan::math::check_less_or_equal;
      check_nonnegative("choose", "n", n);
      check_nonnegative("choose", "k", k);
      if (k > n) return 0;
      const double choices = boost::math::binomial_coefficient<double>(n, k);
      check_less_or_equal("choose", "n choose k", choices,
                          std::numeric_limits<int>::max());
      return std::floor(choices);
    }

  }
}
#endif

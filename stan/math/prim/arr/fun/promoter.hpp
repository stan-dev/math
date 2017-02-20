#ifndef STAN_MATH_PRIM_ARR_FUN_PROMOTER_HPP
#define STAN_MATH_PRIM_ARR_FUN_PROMOTER_HPP

#include <stan/math/prim/scal/fun/promoter.hpp>
#include <vector>

namespace stan {
  namespace math {
    /**
     * Struct to hold static function for element-wise type promotion.
     * This specialization is for vectors of differing types.
     *
     * @tparam F type of input element
     * @tparam T type of output element
     */
    template <typename F, typename T>
    struct promoter<std::vector<F>, std::vector<T> > {
      inline static void promote(const std::vector<F>& u,
                          std::vector<T>& t) {
        t.resize(u.size());
        for (size_t i = 0; i < u.size(); ++i)
          promoter<F, T>::promote(u[i], t[i]);
      }
      inline static std::vector<T>
      promote_to(const std::vector<F>& u) {
        std::vector<T> t;
        promoter<std::vector<F>, std::vector<T> >::promote(u, t);
        return t;
      }
    };

    /**
     * Struct to hold static function for element-wise type promotion.
     * This specialization is for vectors of same types.
     *
     * @tparam F type of input element
     * @tparam T type of output element
     */
    template <typename T>
    struct promoter<std::vector<T>, std::vector<T> > {
      inline static void promote(const std::vector<T>& u,
                          std::vector<T>& t) {
        t = u;
      }
      inline static std::vector<T> promote_to(const std::vector<T>& u) {
        return u;
      }
    };


  }
}

#endif

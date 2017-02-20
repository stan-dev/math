#ifndef STAN_MATH_PRIM_SCAL_FUN_PROMOTER_HPP
#define STAN_MATH_PRIM_SCAL_FUN_PROMOTER_HPP


namespace stan {
  namespace math {
    /**
     * Struct which holds static function for element-wise type promotion.
     *
     * @tparam F type of input element
     * @tparam T type of output element
     */
    template <typename F, typename T>
    struct promoter {
      inline static void promote(const F& u, T& t) {
        t = u;
      }
      inline static T promote_to(const F& u) {
        return u;
      }
    };

    /**
     * Struct which holds static function for element-wise type promotion.
     * This specialization is for scalars of the same type.
     *
     * @tparam T type of input, output element
     */
    template <typename T>
    struct promoter<T, T> {
      inline static void promote(const T& u, T& t) {
        t = u;
      }
      inline static T promote_to(const T& u) {
        return u;
      }
    };

  }
}

#endif

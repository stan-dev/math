#ifndef STAN_MATH_REV_CORE_VARI_VISITOR_HPP
#define STAN_MATH_REV_CORE_VARI_VISITOR_HPP

#include <stan/math/rev/core/vari.hpp>
#include <utility>

namespace stan {
namespace math {

template <size_t N, typename = void>
struct visit_nt_impl;

template <size_t N>
struct visit_nt_impl<N, require_t<bool_constant<(N == 0)>>> {
  template <typename Variant, typename Visitor>
  inline static constexpr auto apply(Variant &&var, Visitor &&vis) noexcept {
        if (N == var.index()) {
        // If this check isnt there the compiler will generate
        // exception code, this stops that
        return std::forward<Visitor>(vis)(
            boost::variant2::get<N>(std::forward<Variant>(var)));
        }
  }
};
template <size_t N>
struct visit_nt_impl<N, require_t<bool_constant<(N > 0)>>> {
  template <typename Variant, typename Visitor>
  inline static constexpr auto apply(Variant &&var, Visitor &&vis) noexcept {
        if (likely(var.index() == N)) {
        return std::forward<Visitor>(vis)(
            *boost::variant2::get_if<vari_value<double>*>(&std::forward<Variant>(var)));
        } else {
        return visit_nt_impl<N - 1>::apply(std::forward<Variant>(var),
                                std::forward<Visitor>(vis));
        }
  }
};

template <typename Variant, typename Visitor>
constexpr auto vari_visitor(Variant&& var, Visitor&& vis) noexcept {
  constexpr auto sizer = boost::variant2::variant_size<Variant>::value;
  return visit_nt_impl<sizer - 1>::apply(std::forward<Variant>(var), std::forward<Visitor>(vis));
}

}
}
#endif

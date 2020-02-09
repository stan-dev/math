#ifndef STAN_MATH_REV_CORE_OP_VARI_HPP
#define STAN_MATH_REV_CORE_OP_VARI_HPP

#include <stan/math/rev/core/vari.hpp>
#include <tuple>

namespace stan {
namespace math {

template <typename... Types>
class op_vari : public vari {
 protected:
    std::tuple<Types...> vi_;
 public:
   template <std::size_t Ind>
   auto& get() {
     return std::get<Ind>(vi_);
   }

  op_vari(double f, Types... args)
      : vari(f), vi_(std::make_tuple(args...)) {}
};

}  // namespace math
}  // namespace stan
#endif

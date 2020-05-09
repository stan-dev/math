#ifndef STAN_MATH_REV_CORE_OP_VARI_HPP
#define STAN_MATH_REV_CORE_OP_VARI_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/prim/meta.hpp>
#include <tuple>

namespace stan {
namespace math {

/**
 * Holds the elements needed in vari operations for the reverse pass and chain
 * call.
 *
 * @tparam Types The types of the operation.
 */
template <typename T, typename... Types>
class op_vari : public vari_value<T> {
 protected:
  std::tuple<Types...> vi_;  // Holds the objects needed in the reverse pass.

 public:
  /**
   * Get an element from the tuple of vari ops. Because of name lookup rules
   *  this function needs to be called as \c this->template get<N>()
   * @tparam Ind The index of the tuple to retrieve.
   * @return the element inside of the tuple at index Ind.
   */
  template <std::size_t Ind>
  auto& get() {
    return std::get<Ind>(vi_);
  }

  /**
   * Return a reference to the tuple holding the vari ops. This is commonly
   *  used in conjunction with \c std::get<N>()
   * @return The tuple holding the vari ops.
   */
  auto& vi() { return vi_; }
  auto& avi() { return std::get<0>(vi_); }
  auto& ad() { return std::get<0>(vi_); }
  auto& bvi() { return std::get<1>(vi_); }
  auto& bd() { return std::get<1>(vi_); }
  auto& cvi() { return std::get<2>(vi_); }
  auto& cd() { return std::get<2>(vi_); }
  /**
   * Constructor for passing in vari and ops objects.
   * @param val Value to initialize the vari to.
   * @param args Ops passed into the tuple and used later in chain method.
   */
  op_vari(T val, Types... args)
      : vari_value<T>(val), vi_(std::make_tuple(args...)) {}
};

}  // namespace math
}  // namespace stan
#endif

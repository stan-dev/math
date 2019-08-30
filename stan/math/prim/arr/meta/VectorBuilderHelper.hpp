#ifndef STAN_MATH_PRIM_ARR_META_VECTORBUILDER_HELPER_HPP
#define STAN_MATH_PRIM_ARR_META_VECTORBUILDER_HELPER_HPP

#include <stan/math/prim/scal/meta/VectorBuilderHelper.hpp>
#include <stdexcept>
#include <vector>
#include <utility>

namespace stan {

/**
 * Template specialization for using a vector
 */
template <typename T1>
class VectorBuilderHelper<T1, true, true> {  // When used and vector
 private:
  std::vector<T1> x_;

 public:
  explicit VectorBuilderHelper(size_t n) : x_(std::move(n)) {}

  typedef std::vector<T1> type;

  auto& operator[](size_t i) { return x_[i]; }
  const auto& operator[](size_t i) const { return x_[i]; }
  inline auto& data() { return x_; }
  inline const auto& data() const { return x_; }
};
}  // namespace stan
#endif

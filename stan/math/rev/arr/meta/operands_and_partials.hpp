#ifndef STAN_MATH_REV_ARR_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_REV_ARR_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <stan/math/prim/scal/meta/likely.hpp>
#include <stan/math/rev/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/arr/meta/length.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {
// Vectorized Univariate
template <>
class ops_partials_edge<double, std::vector<var> > {
 public:
  typedef std::vector<var> Op;
  typedef std::vector<double> partials_t;
  partials_t partials_;                       // For univariate use-cases
  broadcast_array<partials_t> partials_vec_;  // For multivariate
  explicit ops_partials_edge(const Op& op)
      : partials_(partials_t(op.size(), 0)),
        partials_vec_(partials_),
        operands_(op) {}

 private:
  template <typename, typename, typename, typename, typename, typename>
  friend class stan::math::operands_and_partials;
  const Op& operands_;

  void dump_partials(double* partials) {
    for (int i = 0; i < this->partials_.size(); ++i) {
      partials[i] = this->partials_[i];
    }
  }
  void dump_operands(vari** varis) {
    for (size_t i = 0; i < this->operands_.size(); ++i) {
      varis[i] = this->operands_[i].vi_;
    }
  }
  int size() { return this->operands_.size(); }
};
}  // namespace internal
}  // namespace math
}  // namespace stan

#endif

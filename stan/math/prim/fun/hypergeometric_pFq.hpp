#ifndef STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_PFQ_HPP
#define STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

namespace stan {
namespace math {

double hypergeometric_pFq(const Eigen::VectorXd& p, const Eigen::VectorXd& q, double z) {

  return boost::math::hypergeometric_pFq(std::vector<double>(p.data(),p.data()+p.size()),
                                         std::vector<double>(q.data(),q.data()+q.size()),
                                         z);
}
}  // namespace math
}  // namespace stan
#endif

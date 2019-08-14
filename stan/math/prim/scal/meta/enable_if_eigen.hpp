#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_EIGEN_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_EIGEN_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

template <typename T>
using enable_if_eigen = std::enable_if_t<std::is_base_of<Eigen::EigenBase<std::decay_t<T>>, std::decay_t<T>>::value>;

template <typename T>
using enable_if_not_eigen
    = std::enable_if_t<!std::is_base_of<Eigen::EigenBase<std::decay_t<T>>, std::decay_t<T>>::value>;

template <typename... Types>
using enable_if_all_eigen
    = std::enable_if_t<math::conjunction<std::is_base_of<Eigen::EigenBase<std::decay_t<Types>>, std::decay_t<Types>>...>::value>;

template <typename... Types>
using enable_if_any_eigen
    = std::enable_if_t<math::disjunction<std::is_base_of<Eigen::EigenBase<std::decay_t<Types>>, std::decay_t<Types>>...>::value>;

template <typename... Types>
using enable_if_all_not_eigen
    = std::enable_if_t<!math::conjunction<std::is_base_of<Eigen::EigenBase<std::decay_t<Types>>, std::decay_t<Types>>...>::value>;

template <typename... Types>
using enable_if_any_not_eigen
    = std::enable_if_t<!math::disjunction<std::is_base_of<Eigen::EigenBase<std::decay_t<Types>>, std::decay_t<Types>>...>::value>;

}  // namespace stan
#endif

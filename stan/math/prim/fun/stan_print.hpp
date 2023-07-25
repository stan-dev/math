#ifndef STAN_MATH_PRIM_FUN_STAN_PRINT_HPP
#define STAN_MATH_PRIM_FUN_STAN_PRINT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <vector>

namespace stan {
namespace math {
// prints used in generator for print() statements in modeling language

template <typename T, require_not_container_t<T>* = nullptr,
          require_not_tuple_t<T>* = nullptr>
void stan_print(std::ostream* o, const T& x) {
  *o << x;
}

template <typename EigVec, require_eigen_vector_t<EigVec>* = nullptr>
void stan_print(std::ostream* o, const EigVec& x) {
  const auto& x_ref = to_ref(x);

  *o << '[';
  for (int i = 0; i < x_ref.size(); ++i) {
    if (i > 0) {
      *o << ',';
    }
    stan_print(o, x_ref.coeff(i));
  }
  *o << ']';
}

template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_not_eigen_vector_t<EigMat>* = nullptr>
void stan_print(std::ostream* o, const EigMat& x) {
  const auto& x_ref = to_ref(x);

  *o << '[';
  for (int i = 0; i < x_ref.rows(); ++i) {
    if (i > 0) {
      *o << ',';
    }
    *o << '[';
    for (int j = 0; j < x_ref.cols(); ++j) {
      if (j > 0) {
        *o << ',';
      }
      stan_print(o, x_ref.coeff(i, j));
    }
    *o << ']';
  }
  *o << ']';
}

// forward decl to allow the next two overloads to call each other
template <typename T, require_tuple_t<T>* = nullptr>
void stan_print(std::ostream* o, const T& x);

template <typename T, require_std_vector_t<T>* = nullptr>
void stan_print(std::ostream* o, const T& x) {
  *o << '[';
  for (size_t i = 0; i < x.size(); ++i) {
    if (i > 0) {
      *o << ',';
    }
    stan_print(o, x[i]);
  }
  *o << ']';
}

template <typename T, require_tuple_t<T>*>
void stan_print(std::ostream* o, const T& x) {
  *o << '(';
  constexpr auto tuple_size = std::tuple_size<std::decay_t<T>>::value;
  size_t i = 0;
  stan::math::for_each(
      [&i, o](auto&& elt) {
        if (i > 0) {
          *o << ',';
        }
        stan_print(o, elt);
        i++;
      },
      x);
  *o << ')';
}

}  // namespace math
}  // namespace stan
#endif

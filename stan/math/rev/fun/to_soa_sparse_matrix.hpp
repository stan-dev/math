#ifndef STAN_MATH_REV_FUN_TO_SOA_SPARSE_MATRIX_HPP
#define STAN_MATH_REV_FUN_TO_SOA_SPARSE_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
namespace math {

template <int Options = Eigen::ColMajor, typename T1, typename T2, typename T3,
  require_var_t<T1>* = nullptr,
  require_eigen_dense_base_t<value_type_t<T1>>* = nullptr,
  require_all_std_vector_vt<std::is_integral, T2, T3>* = nullptr>
inline auto to_soa_sparse_matrix(int m, int n, T1&& w, T2&& u, T3&& v) {
    arena_t<T2> v_arena(std::forward<T2>(v));
    arena_t<T3> u_arena(std::forward<T3>(u));
    using sparse_arena_mat_t = arena_t<Eigen::SparseMatrix<double, Options>>;
    sparse_arena_mat_t arena_val_x(m, n, w.val().size(),
      u_arena.data(), v_arena.data(), w.vi_->val_.data());
    var_value<Eigen::SparseMatrix<double, Options>> var_x(arena_val_x);
    sparse_arena_mat_t arena_adj_x(m, n, w.adj().size(),
      u_arena.data(), v_arena.data(), w.vi_->adj_.data());
    reverse_pass_callback([arena_adj_x, var_x]() mutable {
      std::cout << "pre varx adj: \n" << var_x.adj() << std::endl;
      std::cout << "pre arena_adj: \n" << arena_adj_x << std::endl;
      arena_adj_x += var_x.adj();
      std::cout << "post varx adj: \n" << var_x.adj() << std::endl;
      std::cout << "post arena_adj: \n" << arena_adj_x << std::endl;
    });
    return var_x;
  }


template <int Options = Eigen::ColMajor, typename T1, typename T2, typename T3,
  require_eigen_dense_base_vt<is_var, T1>* = nullptr,
  require_all_std_vector_vt<std::is_integral, T2, T3>* = nullptr>
inline auto to_soa_sparse_matrix(int m, int n, T1&& w, T2&& u, T3&& v) {
    arena_t<T1> w_arena(std::forward<T1>(w));
    arena_t<T2> v_arena(std::forward<T2>(v));
    arena_t<T3> u_arena(std::forward<T3>(u));
    arena_t<Eigen::SparseMatrix<var, Options>> arena_x(m, n, w_arena.size(),
      u_arena.data(), v_arena.data(), w_arena.data());
    var_value<Eigen::SparseMatrix<double, Options>> var_x(value_of(arena_x));
    // No need to copy adj, but need to backprop
    reverse_pass_callback([arena_x, var_x]() mutable {
      using var_sparse_iterator_t = typename arena_t<Eigen::SparseMatrix<var, Options>>::InnerIterator;
      using dbl_sparse_iterator_t = typename arena_t<Eigen::SparseMatrix<double, Options>>::InnerIterator;
      // arena_x.adj() += var_x.adj() once custom adj() for var sparse matrix
      for (int k = 0; k < arena_x.outerSize(); ++k) {
        var_sparse_iterator_t it_arena_x(arena_x, k);
        dbl_sparse_iterator_t it_var_x(var_x.adj(), k);
        for (; bool(it_arena_x) && bool(it_var_x); ++it_arena_x, ++it_var_x) {
          it_arena_x.valueRef().adj() += it_var_x.valueRef();
        }
      }
    });
    return var_x;
  }

template <int Options = Eigen::ColMajor, typename T1, typename T2, typename T3,
  require_eigen_dense_base_vt<std::is_arithmetic, T1>* = nullptr,
  require_all_std_vector_vt<std::is_integral, T2, T3>* = nullptr>
  inline auto to_soa_sparse_matrix(int m, int n, T1&& w, T2&& u, T3&& v) {
    arena_t<T1> w_arena(std::forward<T1>(w));
    arena_t<T2> v_arena(std::forward<T2>(v));
    arena_t<T3> u_arena(std::forward<T3>(u));
    arena_t<Eigen::SparseMatrix<double, Options>> arena_x(m, n, w_arena.size(),
      u_arena.data(), v_arena.data(), w_arena.data());
    return var_value<Eigen::SparseMatrix<double, Options>>(arena_x);
  }

  template <int Options = Eigen::ColMajor, typename T, require_t<math::conjunction<is_var<T>, is_eigen_sparse_base<value_type_t<T>>>>* = nullptr>
  inline auto to_soa_sparse_matrix(T&& x) {
    return std::forward<T>(x);
  }

}
}

#endif


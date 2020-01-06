#ifndef STAN_MATH_REV_SCAL_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_REDUCE_SUM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>

#include <tbb/task_arena.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

#include <iostream>
#include <iterator>
#include <vector>

namespace stan {
namespace math {

template <typename T>
struct value_type {
  using type = double;
};

template <>
struct value_type<int> {
  using type = int;
};

template <template <typename> class T_struct, typename T_value>
struct value_type<T_struct<T_value>> {
  using type = T_struct<typename value_type<T_value>::type>;
};

template <typename T_value>
struct value_type<std::vector<T_value>> {
  using type = std::vector<typename value_type<T_value>::type>;
};

template <typename T_value, int R, int C>
struct value_type<Eigen::Matrix<T_value, R, C>> {
  using type = Eigen::Matrix<typename value_type<T_value>::type, R, C>;
};

// adjoint_of grabs the adjoint of some data structure

template <typename T>
typename value_type<T>::type adjoint_of(const T x) {
  return typename value_type<T>::type(0);
}

template <typename T>
typename value_type<std::vector<T>>::type adjoint_of(const std::vector<T>& x) {
  const size_t size = x.size();
  typename value_type<std::vector<T>>::type result(size);
  for (size_t i = 0; i != size; ++i) {
    result[i] = adjoint_of(x[i]);
  }
  return result;
}

double adjoint_of(const var& x) { return x.vi_->adj_; }

template <typename T, int R, int C>
inline Eigen::Matrix<typename value_type<T>::type, R, C> adjoint_of(
    const Eigen::Matrix<T, R, C>& M) {
  Eigen::Matrix<typename value_type<T>::type, R, C> Md(M.rows(), M.cols());
  for (int j = 0; j < M.cols(); j++) {
    for (int i = 0; i < M.rows(); i++) {
      Md(i, j) = adjoint_of(M(i, j));
    }
  }
  return Md;
}

template <int R, int C>
inline const Eigen::Matrix<double, R, C>& adjoint_of(
    const Eigen::Matrix<double, R, C>& x) {
  return Eigen::Matrix<double, R, C>::Zero(x.rows(), x.cols());
}

template <int R, int C>
inline const Eigen::Matrix<int, R, C>& adjoint_of(
    const Eigen::Matrix<int, R, C>& x) {
  return Eigen::Matrix<int, R, C>::Zero(x.rows(), x.cols());
}

// as_value is similar to value_of... with the difference that it
// keeps nested container structures which is not the case for
// value_of

template <typename T>
typename value_type<T>::type as_value(const T x) {
  return typename value_type<T>::type(x);
}

template <typename T>
typename value_type<std::vector<T>>::type as_value(const std::vector<T>& x) {
  const size_t size = x.size();
  typename value_type<std::vector<T>>::type result(size);
  for (size_t i = 0; i != size; ++i) {
    result[i] = as_value(x[i]);
  }
  return result;
}

double as_value(const var& x) { return x.vi_->val_; }

template <typename T, int R, int C>
inline Eigen::Matrix<typename value_type<T>::type, R, C> as_value(
    const Eigen::Matrix<T, R, C>& M) {
  Eigen::Matrix<typename value_type<T>::type, R, C> Md(M.rows(), M.cols());
  for (int j = 0; j < M.cols(); j++) {
    for (int i = 0; i < M.rows(); i++) {
      Md(i, j) = as_value(M(i, j));
    }
  }
  return Md;
}

template <int R, int C>
inline const Eigen::Matrix<double, R, C>& as_value(
    const Eigen::Matrix<double, R, C>& x) {
  return x.val();
}

template <int R, int C>
inline const Eigen::Matrix<int, R, C>& as_value(
    const Eigen::Matrix<int, R, C>& x) {
  return x.val();
}

namespace internal {

template <typename Op>
void add_adjoints(ops_partials_edge<double, std::vector<Op>>& edge,
                  const std::vector<typename value_type<Op>::type>& adjoint,
                  std::size_t offset = 0) {
  for (std::size_t i = 0; i != adjoint.size(); ++i)
    edge.partials_[offset + i] += adjoint[i];
}

template <typename Op>
void sum_adjoints(ops_partials_edge<double, std::vector<Op>>& edge1,
                  const ops_partials_edge<double, std::vector<Op>>& edge2) {
  for (std::size_t i = 0; i != edge1.partials_.size(); ++i)
    edge1.partials_[i] += edge2.partials_[i];
}

template <typename Op>
void add_adjoints(
    ops_partials_edge<double, std::vector<std::vector<Op>>>& edge,
    const std::vector<std::vector<typename value_type<Op>::type>>& adjoint,
    std::size_t offset = 0) {
  for (std::size_t i = 0; i != adjoint.size(); ++i)
    for (std::size_t j = 0; j != adjoint[i].size(); ++j)
      edge.partials_vec_[offset + i][j] += adjoint[i][j];
}

template <typename Op>
void sum_adjoints(
    ops_partials_edge<double, std::vector<std::vector<Op>>>& edge1,
    const ops_partials_edge<double, std::vector<std::vector<Op>>>& edge2) {
  for (std::size_t i = 0; i != edge1.partials_.size(); ++i)
    for (std::size_t j = 0; j != edge1.partials_vec_[i].size(); ++j)
      edge1.partials_vec_[i][j] += edge2.partials_vec_[i][j];
}

template <typename Op, int R, int C>
void add_adjoints(
    ops_partials_edge<double, std::vector<Eigen::Matrix<Op, R, C>>>& edge,
    const std::vector<typename value_type<Eigen::Matrix<Op, R, C>>::type>&
        adjoint,
    std::size_t offset = 0) {
  for (std::size_t i = 0; i != adjoint.size(); ++i)
    for (std::size_t j = 0; j != adjoint[i].size(); ++j)
      edge.partials_vec_[offset + i](j) += adjoint[i](j);
}

template <typename Op, int R, int C>
void sum_adjoints(
    ops_partials_edge<double, std::vector<Eigen::Matrix<Op, R, C>>>& edge1,
    const ops_partials_edge<double, std::vector<Eigen::Matrix<Op, R, C>>>&
        edge2) {
  for (std::size_t i = 0; i != edge1.partials_vec_.size(); ++i)
    for (std::size_t j = 0; j != edge1.partials_vec_[i].size(); ++j)
      edge1.partials_vec_[i](j) += edge2.partials_vec_[i](j);
}

template <>
void add_adjoints(ops_partials_edge<double, std::vector<int>>& edge,
                  const std::vector<typename value_type<int>::type>& adjoint,
                  std::size_t offset) {}

template <>
void sum_adjoints(ops_partials_edge<double, std::vector<int>>& edge1,
                  const ops_partials_edge<double, std::vector<int>>& edge2) {}

template <>
void add_adjoints(
    ops_partials_edge<double, std::vector<std::vector<int>>>& edge,
    const std::vector<std::vector<typename value_type<int>::type>>& adjoint,
    std::size_t offset) {}

template <>
void sum_adjoints(
    ops_partials_edge<double, std::vector<std::vector<int>>>& edge1,
    const ops_partials_edge<double, std::vector<std::vector<int>>>& edge2) {}

template <int R, int C>
void add_adjoints(
    ops_partials_edge<double, std::vector<Eigen::Matrix<int, R, C>>>& edge,
    const std::vector<typename value_type<Eigen::Matrix<int, R, C>>::type>&
        adjoint,
    std::size_t offset) {}

template <int R, int C>
void sum_adjoints(
    ops_partials_edge<double, std::vector<Eigen::Matrix<int, R, C>>>& edge1,
    const ops_partials_edge<double, std::vector<Eigen::Matrix<int, R, C>>>&
        edge2) {}

template <class ReduceFunction, class M, class T, class Arg1, class Arg2,
          class Arg3, class Arg4>
struct reduce_sum_impl<ReduceFunction, M, T, Arg1, Arg2, Arg3, Arg4, var> {
  using vmapped_t = std::vector<M>;
  using arg1_t = std::vector<Arg1>;
  using arg2_t = std::vector<Arg2>;
  using arg3_t = std::vector<Arg3>;
  using arg4_t = std::vector<Arg4>;

  using vmapped_value_t = std::vector<typename value_type<M>::type>;
  using arg1_value_t = std::vector<typename value_type<Arg1>::type>;
  using arg2_value_t = std::vector<typename value_type<Arg2>::type>;
  using arg3_value_t = std::vector<typename value_type<Arg3>::type>;
  using arg4_value_t = std::vector<typename value_type<Arg4>::type>;

  using ops_partials_t
      = operands_and_partials<vmapped_t, std::vector<Arg1>, std::vector<Arg2>,
                              std::vector<Arg3>, std::vector<Arg4>>;

  struct recursive_reducer {
    const vmapped_t& vmapped_;
    const arg1_t& arg1_;
    const arg1_value_t& arg1_value_;
    const arg2_t& arg2_;
    const arg2_value_t& arg2_value_;
    const arg3_t& arg3_;
    const arg3_value_t& arg3_value_;
    const arg4_t& arg4_;
    const arg4_value_t& arg4_value_;

    ops_partials_t& terms_partials_mapped_;
    ops_partials_t terms_partials_args_;
    double terms_sum_;

    recursive_reducer(const vmapped_t& vmapped, const T& init,
                      const arg1_t& arg1, const arg1_value_t& arg1_value,
                      const arg2_t& arg2, const arg2_value_t& arg2_value,
                      const arg3_t& arg3, const arg3_value_t& arg3_value,
                      const arg4_t& arg4, const arg4_value_t& arg4_value,
                      ops_partials_t& terms_partials_mapped)
        : vmapped_(vmapped),
          arg1_(arg1),
          arg1_value_(arg1_value),
          arg2_(arg2),
          arg2_value_(arg2_value),
          arg3_(arg3),
          arg3_value_(arg3_value),
          arg4_(arg4),
          arg4_value_(arg4_value),
          terms_partials_mapped_(terms_partials_mapped),
          terms_partials_args_(vmapped, arg1, arg2, arg3, arg4),
          terms_sum_(as_value(init)) {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : vmapped_(other.vmapped_),
          arg1_(other.arg1_),
          arg1_value_(other.arg1_value_),
          arg2_(other.arg2_),
          arg2_value_(other.arg2_value_),
          arg3_(other.arg3_),
          arg3_value_(other.arg3_value_),
          arg4_(other.arg4_),
          arg4_value_(other.arg4_value_),
          terms_partials_mapped_(other.terms_partials_mapped_),
          terms_partials_args_(vmapped_, arg1_, arg2_, arg3_, arg4_),
          terms_sum_(0.0) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty())
        return;

      try {
        start_nested();

        // create a deep copy of all var's so that these are not
        // linked to any outer AD tree

        vmapped_t local_sub_slice;

        local_sub_slice.reserve(r.size());
        for (std::size_t i = r.begin(); i != r.end(); ++i)
          local_sub_slice.emplace_back(as_value(vmapped_[i]));

        const arg1_t local_arg1(arg1_value_.begin(), arg1_value_.end());
        const arg2_t local_arg2(arg2_value_.begin(), arg2_value_.end());
        const arg3_t local_arg3(arg3_value_.begin(), arg3_value_.end());
        const arg4_t local_arg4(arg4_value_.begin(), arg4_value_.end());

        T sub_sum_v
            = ReduceFunction()(r.begin(), r.end() - 1, local_sub_slice,
                               local_arg1, local_arg2, local_arg3, local_arg4);

        sub_sum_v.grad();

        terms_sum_ += sub_sum_v.val();

        if (!is_constant_all<M>::value) {
          vmapped_value_t sub_slice_adjoint = adjoint_of(local_sub_slice);
          add_adjoints(terms_partials_mapped_.edge1_, sub_slice_adjoint,
                       r.begin());
        }

        if (!is_constant_all<Arg1>::value) {
          arg1_value_t arg1_adjoint = adjoint_of(local_arg1);
          add_adjoints(terms_partials_args_.edge2_, arg1_adjoint);
        }
        if (!is_constant_all<Arg2>::value) {
          arg2_value_t arg2_adjoint = adjoint_of(local_arg2);
          add_adjoints(terms_partials_args_.edge3_, arg2_adjoint);
        }
        if (!is_constant_all<Arg3>::value) {
          arg3_value_t arg3_adjoint = adjoint_of(local_arg3);
          add_adjoints(terms_partials_args_.edge4_, arg3_adjoint);
        }
        if (!is_constant_all<Arg4>::value) {
          arg4_value_t arg4_adjoint = adjoint_of(local_arg4);
          add_adjoints(terms_partials_args_.edge5_, arg4_adjoint);
        }

      } catch (const std::exception& e) {
        recover_memory_nested();
        throw;
      }
      recover_memory_nested();
    }

    void join(const recursive_reducer& child) {
      terms_sum_ += child.terms_sum_;

      if (!is_constant_all<Arg1>::value) {
        sum_adjoints(terms_partials_args_.edge2_,
                     child.terms_partials_args_.edge2_);
      }
      if (!is_constant_all<Arg2>::value) {
        sum_adjoints(terms_partials_args_.edge3_,
                     child.terms_partials_args_.edge3_);
      }
      if (!is_constant_all<Arg3>::value) {
        sum_adjoints(terms_partials_args_.edge4_,
                     child.terms_partials_args_.edge4_);
      }
      if (!is_constant_all<Arg4>::value) {
        sum_adjoints(terms_partials_args_.edge5_,
                     child.terms_partials_args_.edge5_);
      }
    }
  };

  T operator()(const vmapped_t& vmapped, T init, std::size_t grainsize,
               const arg1_t& arg1, const arg2_t& arg2, const arg3_t& arg3,
               const arg4_t& arg4) const {
    const std::size_t num_jobs = vmapped.size();
    const arg1_value_t arg1_value = as_value(arg1);
    const arg2_value_t arg2_value = as_value(arg2);
    const arg3_value_t arg3_value = as_value(arg3);
    const arg4_value_t arg4_value = as_value(arg4);

    ops_partials_t ops(vmapped, arg1, arg2, arg3, arg4);

    recursive_reducer worker(vmapped, init, arg1, arg1_value, arg2, arg2_value,
                             arg3, arg3_value, arg4, arg4_value, ops);

    tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker);

    if (!is_constant_all<Arg1>::value) {
      ops.edge2_.partials_vec_
          = worker.terms_partials_args_.edge2_.partials_vec_;
    }
    if (!is_constant_all<Arg2>::value) {
      ops.edge3_.partials_vec_
          = worker.terms_partials_args_.edge3_.partials_vec_;
    }
    if (!is_constant_all<Arg3>::value) {
      ops.edge4_.partials_vec_
          = worker.terms_partials_args_.edge4_.partials_vec_;
    }
    if (!is_constant_all<Arg4>::value) {
      ops.edge5_.partials_vec_
          = worker.terms_partials_args_.edge5_.partials_vec_;
    }

    return ops.build(worker.terms_sum_);
  }
};
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif

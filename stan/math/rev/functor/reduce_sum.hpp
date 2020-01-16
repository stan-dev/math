#ifndef STAN_MATH_REV_SCAL_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_REDUCE_SUM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/fun/typedefs.hpp>

#include <tbb/task_arena.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

#include <iostream>
#include <iterator>
#include <vector>

namespace stan {
namespace math {
namespace internal {

template <typename ReduceFunction, typename ReturnType, typename M,
          typename... Args>
struct reduce_sum_impl<ReduceFunction, require_var_t<ReturnType>, ReturnType, M,
                       Args...> {
  template <typename T>
  static T& deep_copy(T& arg) {
    return arg;
  }

  static var deep_copy(const var& arg) { return var(arg.val()); }

  static std::vector<var> deep_copy(const std::vector<var>& arg) {
    std::vector<var> copy(arg.size());
    for (size_t i = 0; i < arg.size(); ++i) {
      copy[i] = arg[i].val();
    }
    return copy;
  }

  template <int RowType, int ColType>
  static Eigen::Matrix<var, RowType, ColType> deep_copy(
      const Eigen::Matrix<var, RowType, ColType>& arg) {
    Eigen::Matrix<var, RowType, ColType> copy(arg.size());
    for (size_t i = 0; i < arg.size(); ++i) {
      copy(i) = arg(i).val();
    }
    return copy;
  }

  // TODO(Steve): Move this somewhere smarter
  // Fails to compile if type T does not have member operator(Integral)
  template <typename T>
  using operator_paren_access_t = decltype(std::declval<T>()(int{}));

  template <typename... Pargs>
  static double* accumulate_adjoints(double* dest, const var& x,
                                     const Pargs&... args) {
    *dest += x.adj();
    return accumulate_adjoints(dest + 1, args...);
  }

  // Works with anything that has operator[Integral] defined
  template <typename... Pargs, typename Vec,
            require_vector_like_vt<is_var, Vec>...>
  static double* accumulate_adjoints(double* dest, const Vec& x,
                                     const Pargs&... args) {
    for (size_t i = 0; i < x.size(); ++i) {
      dest[i] += x[i].adj();
    }
    return accumulate_adjoints(dest + x.size(), args...);
  }

  // Works on anything with a operator()
  template <typename... Pargs, typename Mat,
            require_t<is_detected<Mat, operator_paren_access_t>>...,
            require_t<is_var<value_type_t<Mat>>>...>
  static double* accumulate_adjoints(double* dest, const Mat& x,
                                     const Pargs&... args) {
    Eigen::Map<matrix_vi>(dest, x.size()) += x.adj();
    return accumulate_adjoints(dest + x.size(), args...);
  }

  // Anything with a scalar type of Arithmetic gets tossed
  template <typename Arith, require_arithmetic_t<scalar_type_t<Arith>>...,
            typename... Pargs>
  static double* accumulate_adjoints(double* dest, Arith&& x,
                                     const Pargs&... args) {
    return accumulate_adjoints(dest, args...);
  }

  static double* accumulate_adjoints(double* x) { return x; }

  struct recursive_reducer {
    size_t num_shared_terms_;
    const std::vector<M>& vmapped_;
    std::tuple<const Args&...> args_tuple_;
    size_t tuple_size_ = sizeof...(Args);

    double sum_;
    Eigen::VectorXd args_adjoints_;

    recursive_reducer(size_t num_shared_terms, const std::vector<M>& vmapped,
                      const Args&... args)
        : num_shared_terms_(num_shared_terms),
          vmapped_(vmapped),
          args_tuple_(args...),
          sum_(0.0),
          args_adjoints_(Eigen::VectorXd::Zero(num_shared_terms_)) {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : num_shared_terms_(other.num_shared_terms_),
          vmapped_(other.vmapped_),
          args_tuple_(other.args_tuple_),
          sum_(0.0),
          args_adjoints_(Eigen::VectorXd::Zero(num_shared_terms_)) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty())
        return;

      auto start = vmapped_.begin();
      std::advance(start, r.begin());
      auto end = vmapped_.begin();
      std::advance(end, r.end());

      try {
        start_nested();

        // create a deep copy of all var's so that these are not
        // linked to any outer AD tree
        const std::vector<M> sub_slice(start, end);

        auto args_tuple_local_copy = apply(
            [&](auto&&... args) { return std::make_tuple(deep_copy(args)...); },
            args_tuple_);

        var sub_sum_v = apply(
            [&r, &sub_slice](auto&&... args) {
              return ReduceFunction()(r.begin(), r.end() - 1, sub_slice,
                                      args...);
            },
            args_tuple_local_copy);

        sub_sum_v.grad();

        sum_ += sub_sum_v.val();

        // This should accumulate the adjoints from args_tuple_local_copy into
        //  the memory of args_adjoints_
        apply(
            [&](auto&&... args) {
              return accumulate_adjoints(args_adjoints_.data(), args...);
            },
            args_tuple_local_copy);
      } catch (const std::exception& e) {
        recover_memory_nested();
        throw;
      }
      recover_memory_nested();
    }

    void join(const recursive_reducer& rhs) {
      sum_ += rhs.sum_;
      args_adjoints_ += rhs.args_adjoints_;
    }
  };

  // Fails to compile if type T does not have callable member size()
  template <typename T>
  using member_size_t = decltype(std::declval<T>().size());

  // TODO(Steve): add requires for generic containers
  template <typename Container,
            require_t<is_detected<Container, member_size_t>>...,
            require_t<is_var<value_type_t<Container>>>..., typename... Pargs>
  size_t count_var_impl(size_t count, const Container& x,
                        const Pargs&... args) const {
    return count_var_impl(count + x.size(), args...);
  }

  // TODO(Steve): add this back if you want it cause Ben commented it out cause
  //  it was causing ambiguities
  /*template <typename Container,
            require_t<is_detected<Container, member_size_t>>...,
            require_t<std::is_arithmetic<value_type_t<Container>>>...,
            typename... Pargs>
  size_t count_var_impl(size_t count, const Container& x,
                        const Pargs&... args) const {
    return count_var_impl(count, args...);
    }*/

  template <typename... Pargs>
  size_t count_var_impl(size_t count, const var& x,
                        const Pargs&... args) const {
    return count_var_impl(count + 1, args...);
  }

  template <typename... Pargs, typename Arith,
            require_arithmetic_t<scalar_type_t<Arith>>...>
  size_t count_var_impl(size_t count, Arith& x, const Pargs&... args) const {
    return count_var_impl(count, args...);
  }

  size_t count_var_impl(size_t count) const { return count; }

  /**
   * Count the number of scalars of type T in the input argument list
   *
   * @tparam Pargs Types of input arguments
   * @return Number of scalars of type T in input
   */
  template <typename... Pargs>
  size_t count_var(const Pargs&... args) const {
    return count_var_impl(0, args...);
  }

  template <typename... Pargs>
  void save_varis(vari** dest, const var& x, const Pargs&... args) const {
    *dest = x.vi_;
    save_varis(dest + 1, args...);
  }

  template <typename... Pargs, typename Vec,
            require_vector_like_vt<is_var, Vec>...>
  void save_varis(vari** dest, const Vec& x, const Pargs&... args) const {
    for (size_t i = 0; i < x.size(); ++i) {
      dest[i] = x[i].vi_;
    }
    save_varis(dest + x.size(), args...);
  }

  template <typename... Pargs, typename Mat,
            require_t<is_detected<Mat, operator_paren_access_t>>...,
            require_t<is_var<value_type_t<Mat>>>...>
  void save_varis(vari** dest, const Mat& x, const Pargs&... args) const {
    for (size_t i = 0; i < x.size(); ++i) {
      dest[i] = x(i).vi_;
    }
    save_varis(dest + x.size(), args...);
  }

  template <typename R, require_arithmetic_t<scalar_type_t<R>>...,
            typename... Pargs>
  void save_varis(vari** dest, const R& x, const Pargs&... args) const {
    save_varis(dest, args...);
  }

  void save_varis(vari**) const {}

  var operator()(const std::vector<M>& vmapped, std::size_t grainsize,
                 const Args&... args) const {
    const std::size_t num_jobs = vmapped.size();

    if (num_jobs == 0)
      return var(0.0);

    const std::size_t num_sliced_terms = count_var(vmapped);
    const std::size_t num_shared_terms = count_var(args...);

    vari** varis = ChainableStack::instance_->memalloc_.alloc_array<vari*>(
        num_sliced_terms + num_shared_terms);
    double* partials = ChainableStack::instance_->memalloc_.alloc_array<double>(
        num_sliced_terms + num_shared_terms);

    double sum = 0;

    try {
      start_nested();

      auto vmapped_copy = deep_copy(vmapped);

      recursive_reducer worker(num_shared_terms, vmapped_copy, args...);

#ifdef STAN_DETERMINISTIC
      tbb::static_partitioner partitioner;
      tbb::parallel_deterministic_reduce(
          tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker,
          partitioner);
#else
      tbb::parallel_reduce(
          tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker);
#endif

      save_varis(varis, vmapped);
      save_varis(varis + num_sliced_terms, args...);

      for (size_t i = 0; i < num_sliced_terms; ++i)
        partials[i] = 0.0;
      accumulate_adjoints(partials, vmapped_copy);

      for (size_t i = 0; i < num_shared_terms; ++i) {
        partials[num_sliced_terms + i] = worker.args_adjoints_(i);
      }

      sum = worker.sum_;

    } catch (const std::exception& e) {
      recover_memory_nested();
      throw;
    }
    recover_memory_nested();

    return var(new precomputed_gradients_vari(
        sum, num_sliced_terms + num_shared_terms, varis, partials));
  }
};
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif

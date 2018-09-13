#ifndef TEST_UNIT_MATH_REV_MAT_UTIL_TEST_AUTODIFF_HPP
#define TEST_UNIT_MATH_REV_MAT_UTIL_TEST_AUTODIFF_HPP

#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/scal/functor/call_all_argument_combos.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace stan {
namespace math {
namespace test {

/**
 * If x1 and x2 have absolute value less than 1, return true if abs(x1 - x2) <
 * absolute_tolerance. Otherwise, return true if abs(x1 - x2) / min(abs(x1),
 * abs(x2)) < relative_tolerance.
 *
 * In any other case return false.
 *
 * @param x1 First scalar to compare
 * @param x2 Second scalar to compare
 * @param relative_tolerance Relative tolerance
 * @param absolute_tolerance Absolute tolerance
 */
bool is_near(double x1, double x2, double relative_tolerance = 1e-9,
             double absolute_tolerance = 1e-9) {
  if (is_nan(x1) || is_nan(x2)) {
    return is_nan(x1) && is_nan(x2);
  } else if (is_inf(x1) || is_inf(x2)) {
    return x1 == x2;
  } else {
    if (x1 < 1.0 && x2 < 1.0) {
      return abs(x1 - x2) < absolute_tolerance;
    } else {
      return abs(x1 - x2) / std::min(abs(x1), abs(x2)) < relative_tolerance;
    }
  }
}

/**
 * Return an easy to read string describing the type of x
 *
 * @tparam T Type of argument
 * @param x argument
 * @return String formatted type
 */
template <typename T>
std::string type_as_string(const T& x) {
  return typeid(x).name();
}

template <>
std::string type_as_string<var>(const var& x) {
  return "var";
}

template <>
std::string type_as_string<double>(const double& x) {
  return "double";
}

template <>
std::string type_as_string<int>(const int& x) {
  return "int";
}

template <typename T, int R, int C>
std::string type_as_string(const Eigen::Matrix<T, R, C>&) {
  std::stringstream s;
  s << "Eigen::Matrix<" << type_as_string(T()) << ", " << R << ", "
    << "C>";
  return s.str();
}

template <typename T>
std::string type_as_string(const std::vector<T>&) {
  std::stringstream s;
  s << "std::vector<" << type_as_string(T()) << ">";
  return s.str();
}

/**
 * Checks that f called with args and f called with args cast to
 * non-autodiff types either both throw exceptions or both do not.
 * This does not check that exceptions are of the same type
 *
 * @tparam F Type of functor f
 * @tparam Targs Types of arguments to pass to f
 * @param f Functor to test
 * @param args Arguments to call f with
 * @throw runtime_error if exceptions do not match
 */
template <typename F, typename... Targs>
void test_exceptions(F f, const Targs&... args) {
  auto fd = [&f](auto... args) { return f(value_of(args)...); };
  std::stringstream s;

  try {
    fd(args...);
    try {
      f(args...);
    } catch (const std::exception& e2) {
      s << "prim version of function threw no exception but var version "
           "threw: \""
        << e2.what() << "\"";
    } catch (...) {
      s << "prim version of function threw no exception but var version threw "
           "unknown exception";
    }
  } catch (const std::exception& e) {
    try {
      f(args...);
      s << "prim version of function threw with message: \"" << e.what()
        << "\" but var version threw nothing";
    } catch (...) {
    }
  } catch (...) {
    try {
      f(args...);
      s << "prim version of function threw unknown exception but var version "
           "threw nothing";
    } catch (...) {
    }
  }

  std::string error_msg = s.str();
  if (error_msg.size() > 0)
    throw std::runtime_error(error_msg);
}

/**
 * Checks that the values of the output of f called with args called with
 * args cast to non-autodiff types are the same.
 *
 * @tparam F Type of functor f
 * @tparam Targs Types of arguments to pass to f
 * @param f Functor to test
 * @param relative_tolerance Relative tolerance for comparisons
 * @param absolute_tolerance Absolute tolerance for comparisons
 * @param args Arguments to call f with
 * @throw runtime_error if values do not match
 */
template <typename F, typename... Targs>
void test_values(F f, double relative_tolerance, double absolute_tolerance,
                 const Targs&... args) {
  auto fd = [&f](auto... args) { return f(value_of(args)...); };
  auto input = make_variable_adapter<var>(args...);
  auto output = make_variable_adapter<var>(apply(f, input.get_args()));

  if (output.size() > 0) {
    auto output_double
        = make_variable_adapter<double>(apply(fd, input.get_args()));

    if (output_double.size() != output.size()) {
      std::cout << "The non-autodiff output has " << output_double.size()
                << " elements and the reverse mode output has " << output.size()
                << " elements" << std::endl;
    }

    for (size_t i = 0; i < output.size(); ++i) {
      if (!is_near(output_double(i), value_of(output(i)), relative_tolerance,
                   absolute_tolerance)) {
        std::stringstream s;
        s << "The " << i << "th output element from the non-autodiff call ("
          << output_double(i)
          << ") does not match the value from the reverse mode call ("
          << value_of(output(i)) << ")" << std::endl;
        throw std::runtime_error(s.str());
      }
    }
  }
}

/**
 * Checks that the Jacobian of f computed with autodiff
 * is the same (within tolerance) of the Jacobian of f computed
 * with finite differences
 *
 * @tparam F Type of functor f
 * @tparam Targs Types of arguments to pass to f
 * @param f Functor to test
 * @param dx Finite difference delta
 * @param relative_tolerance Relative tolerance for comparisons
 * @param absolute_tolerance Absolute tolerance for comparisons
 * @param args Arguments to call f with
 * @throw runtime_error if gradients do not match to within tolerance
 */
template <typename F, typename... Targs>
void test_gradients(F f, double dx, double relative_tolerance,
                    double absolute_tolerance, const Targs&... args) {
  start_nested();

  try {
    auto fd = [&f](auto... args) { return f(value_of(args)...); };
    auto input = make_variable_adapter<var>(args...);
    auto output = make_variable_adapter<var>(apply(f, input.get_args()));

    Eigen::MatrixXd fd_jac(output.size(), input.size());
    for (size_t j = 0; j < input.size(); ++j) {
      auto fd_input = input;
      fd_input(j) += dx;
      auto output1
          = make_variable_adapter<double>(apply(fd, fd_input.get_args()));
      fd_input(j) -= 2 * dx;
      auto output2
          = make_variable_adapter<double>(apply(fd, fd_input.get_args()));
      for (size_t i = 0; i < output.size(); ++i) {
        fd_jac(i, j) = (output1(i) - output2(i)) / (2.0 * dx);
      }
    }

    Eigen::MatrixXd var_jac(output.size(), input.size());
    for (size_t i = 0; i < output.size(); ++i) {
      // This needs to be set_zero_all_adjoints and not set_zero_adjoints_nested
      // because the input vars are created outside this nesting level
      set_zero_all_adjoints();

      output(i).grad();

      for (size_t j = 0; j < input.size(); ++j) {
        var_jac(i, j) = input(j).adj();
      }
    }

    for (int i = 0; i < fd_jac.rows(); ++i) {
      for (int j = 0; j < fd_jac.cols(); ++j) {
        if (!is_near(fd_jac(i, j), var_jac(i, j), relative_tolerance,
                     absolute_tolerance)) {
          std::stringstream s;
          s << "The (" << i << ", " << j
            << ")th Jacobian element of the finite difference approximation ("
            << fd_jac(i, j)
            << ") does not match the one from the reverse mode calculation ("
            << var_jac(i, j) << ")" << std::endl;
          throw std::runtime_error(s.str());
        }
      }
    }
  } catch (...) {
    recover_memory_nested();

    throw;
  }

  recover_memory_nested();
}

/**
 * Test that the functor f called with the arguments args (which may contain
 * vars)
 *
 * 1. Throws exceptions at the same time as the matching non-autodiff call
 * 2. Produces the same output as the matching non-autodiff call
 * 3. Can produce the same Jacobian (to a tolerance) as one computed via finite
 * differences with the matching non-autodiff f
 *
 * If nothing in args is a var this testing is skipped.
 *
 * @tparam F Type of functor f
 * @tparam Targs Types of arguments to pass to f
 * @param f Functor to test
 * @param dx Finite difference delta
 * @param relative_tolerance Relative tolerance for comparisons
 * @param absolute_tolerance Absolute tolerance for comparisons
 * @param args Arguments to call f with
 * @throw runtime_error if exceptions do not match
 * @throw runtime_error if values do not match
 * @throw runtime_error if gradients do not match
 */
template <typename F, typename... Targs>
void test_var(F f, double dx, double relative_tolerance,
              double absolute_tolerance, const Targs&... args) {
  auto input = make_variable_adapter<var>(args...);

  if (input.size() == 0) {
    return;
  }

  try {
    test_exceptions(f, args...);
    test_values(f, relative_tolerance, absolute_tolerance, args...);
    test_gradients(f, dx, relative_tolerance, absolute_tolerance, args...);
  } catch (std::exception& e) {
    std::stringstream s;
    std::vector<std::string> argument_types = {{type_as_string(args)...}};
    s << "Error testing f(";
    for (size_t i = 0; i < argument_types.size(); ++i) {
      s << argument_types[i];
      if (i < argument_types.size() - 1) {
        s << ", ";
      }
    }
    s << ")" << std::endl << e.what();
    throw std::runtime_error(s.str());
  }
}

/**
 * Test reverse mode autodiff on the function f given the non-autodiff
 * arguments args.
 *
 * args can be of type int, double, std::vector<int>, std::vector<double>, or
 * Eigen::Matrix<double, R, C> (where R and C are valid Eigen Row/Column types)
 *
 * An argument can be promoted to an autodiff type if it has scalar doubles (so
 * double, std::vector<double>, and Eigen::Matrix<double, R, C>).
 *
 * test_autodiff calls test_var for every combination of promoted and
 * non-promoted arguments
 *
 * @tparam F Type of Functor f
 * @tparam Targs Types of arguments of f
 * @param f Functor
 * @param relative_tolerance Relative tolerance for comparisons
 * @param absolute_tolerance Absolute tolerance for comparisons
 * @param args Arguments for f
 */
template <typename F, typename... Targs>
void test_autodiff(F f, double relative_tolerance, double absolute_tolerance,
                   const Targs&... args) {
  start_nested();

  try {
    call_all_argument_combos(
        [&](auto... args) {
          test_var(f, 1e-4, relative_tolerance, absolute_tolerance, args...);
          return 0;
        },
        std::tuple_cat(std::make_tuple(args),
                       promote_double_to<var>(std::make_tuple(args)))...);
  } catch (...) {
    recover_memory_nested();

    throw;
  }

  recover_memory_nested();
}

}  // namespace test
}  // namespace math
}  // namespace stan

#endif

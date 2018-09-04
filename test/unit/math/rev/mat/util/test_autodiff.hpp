#ifndef TEST_UNIT_MATH_REV_MAT_UTIL_TEST_AUTODIFF_HPP
#define TEST_UNIT_MATH_REV_MAT_UTIL_TEST_AUTODIFF_HPP

#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/scal/functor/call_all_argument_combos.hpp>
#include <gtest/gtest.h>
#include <stdexcept>

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

class exception_not_thrown : public std::exception {};

/**
 * expect_exception calls f with args checks, to the best of its ability, that
 * the function throws an exception of type E
 *
 * @tparam F Type of functor f
 * @tparam E Type of exception to catch
 * @tparam Targs Types of arguments to pass to f
 * @param f Functor to test
 * @param e Exception that was previously caught
 * @param args Arguments to call f with
 */
template <typename F, typename E, typename... Targs>
void expect_exception(F f, const E& e, const Targs&... args) {
  try {
    f(args...);
    throw exception_not_thrown();
  } catch (const exception_not_thrown&) {
    std::stringstream s;
    s << "prim version of function threw with message: \"" << e.what()
      << "\" but var version threw nothing";
    throw std::runtime_error(s.str());
  } catch (const E&) {
  } catch (const std::exception& e2) {
    std::stringstream s;
    s << "prim version of function threw: \"" << e.what()
      << "\" but var version threw different error: \"" << e2.what() << "\"";
    throw std::runtime_error(s.str());
  } catch (...) {
    std::stringstream s;
    s << "prim version of function threw: \"" << e.what()
      << "\" but var version threw different error";
    throw std::runtime_error(s.str());
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
 * Test that the functor f called with the arguments args (which may contain
 * vars)
 *
 * 1. Throws the same exceptions as f called with args cast to non-autodiff
 * types. WARNING: This check isn't perfect. It is possible that the
 * non-autodiff and autodiff versions of this function throw different
 * exceptions and this check fails. If exception behavior is critical, use
 * different tests
 * 2. Produces the same values as f called with args cast to non-autodiff types
 * 3. Can produce the same Jacobian (to a tolerance) as one computed via finite
 * differences with the version of f called with non-autodiff args
 *
 * If no variable in args is a var or has a var as a scalar type, this testing
 * is skipped.
 *
 * @tparam F Type of functor f
 * @tparam Targs Types of arguments to pass to f
 * @param f Functor to test
 * @param args Arguments to call f with
 * @throw runtime_error if exceptions do not match
 * @throw runtime_error if values do not match
 * @throw runtime_error if gradients do not match
 */
template <typename F, typename... Targs>
void test_var(F f, double dx, const Targs&... args) {
  start_nested();

  auto fd = [&f](auto... args) { return f(value_of(args)...); };
  auto input = make_variable_adapter<var>(args...);

  if (input.size() == 0) {
    return;
  }

  bool no_exceptions = false;

  try {
    try {
      fd(args...);
      throw exception_not_thrown();
    } catch (const std::domain_error& e) {
      expect_exception(f, e, args...);
    } catch (const std::invalid_argument& e) {
      expect_exception(f, e, args...);
    } catch (const std::out_of_range& e) {
      expect_exception(f, e, args...);
    } catch (const std::system_error& e) {
      expect_exception(f, e, args...);
    } catch (const std::runtime_error& e) {
      expect_exception(f, e, args...);
    } catch (const std::logic_error& e) {
      expect_exception(f, e, args...);
    } catch (const exception_not_thrown& e) {
      try {
        f(args...);
        no_exceptions = true;
      } catch (const std::exception& e2) {
        std::stringstream s;
        s << "prim version of function threw no exception but var version "
             "threw: "
          << e2.what();
        throw std::runtime_error(s.str());
      } catch (...) {
        throw std::runtime_error(
            "prim version of function threw no exception but var version threw "
            "unknown exception");
      }
    } catch (const std::exception& e) {
      expect_exception(f, e, args...);
    } catch (...) {
      try {
        f(args...);
        throw exception_not_thrown();
      } catch (const exception_not_thrown&) {
        std::stringstream s;
        s << "prim version of function threw unknown exception but var version "
             "threw nothing";
        throw std::runtime_error(s.str());
      } catch (...) {
      }
    }

    if (!no_exceptions) {
      return;
    }

    auto output = make_variable_adapter<var>(apply(f, input.get_args()));

    if (output.size() == 0) {
      throw std::runtime_error("The function produced no output");
    }

    auto output_double
        = make_variable_adapter<double>(apply(fd, input.get_args()));

    if (output_double.size() != output.size()) {
      std::cout << "The non-autodiff output has " << output_double.size()
                << " elements and the reverse mode output has " << output.size()
                << " elements" << std::endl;
    }

    for (size_t i = 0; i < output.size(); ++i) {
      if (!is_near(output_double(i), value_of(output(i)))) {
        std::stringstream s;
        s << "The " << i << "th output element from the non-autodiff call ("
          << output_double(i)
          << ") does not match the value from the reverse mode call ("
          << value_of(output(i)) << ")" << std::endl;
        throw std::runtime_error(s.str());
      }
    }

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
      output(i).grad();

      for (size_t j = 0; j < input.size(); ++j) {
        var_jac(i, j) = input(j).adj();
      }

      set_zero_all_adjoints();
    }

    for (int i = 0; i < fd_jac.rows(); ++i) {
      for (int j = 0; j < fd_jac.cols(); ++j) {
        if (!is_near(fd_jac(i, j), var_jac(i, j), dx * dx * 1e2)) {
          std::stringstream s;
          s << "The (" << i << ", " << j
            << ")th Jacobian element of the finite difference approximation ("
            << fd_jac(i, j)
            << ") does not match the one from the reverse mode calcluation ("
            << var_jac(i, j) << ")" << std::endl;
          throw std::runtime_error(s.str());
        }
      }
    }
  } catch (std::exception& e) {
    std::stringstream s;
    std::vector<std::string> argument_types = {{type_as_string(args)...}};
    s << "Error testing gradients of f(";
    for (int i = 0; i < argument_types.size(); ++i) {
      s << argument_types[i];
      if (i < argument_types.size() - 1) {
        s << ", ";
      }
    }
    s << ")" << std::endl << e.what();
    throw std::runtime_error(s.str());
  }

  recover_memory_nested();
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
 * @param args Arguments for f
 */
template <typename F, typename... Targs>
void test_autodiff(F f, const Targs&... args) {
  call_all_argument_combos(
      [&f](auto... args) {
        test_var(f, 1e-4, args...);
        return 0;
      },
      std::tuple_cat(std::make_tuple(args),
                     promote_double_to_T<var>(std::make_tuple(args)))...);
}

}  // namespace test
}  // namespace math
}  // namespace stan

#endif

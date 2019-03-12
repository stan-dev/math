#include <stan/math/rev.hpp>
#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/util.hpp>
#include <limits>
#include <string>
#include <vector>

template <typename T_x1, typename T_x2, typename T_sigma, typename T_l>
std::string pull_msg(std::vector<T_x1> x1, std::vector<T_x2> x2, T_sigma sigma,
                     T_l l) {
  std::string message;
  try {
    stan::math::gp_exponential_cov(x1, x2, sigma, l);
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exception";
  }
  return message;
}

template <typename T_x1, typename T_sigma, typename T_l>
std::string pull_msg(std::vector<T_x1> x1, T_sigma sigma, T_l l) {
  std::string message;
  try {
    stan::math::gp_exponential_cov(x1, sigma, l);
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exception";
  }
  return message;
}

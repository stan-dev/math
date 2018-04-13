#include "newfun.hpp"
#include <iostream>

int sassy() {
  std::cout << "One's: " << &stan::math::chainable_stack.memalloc_ << std::endl;
  return 2;
}

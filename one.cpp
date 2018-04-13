#include "newfun.hpp"
#include <iostream>

int sassy() {
  std::cout << "One's: " << &stan::math::ADStack::instance.memalloc_ << std::endl;
  return 2;
}

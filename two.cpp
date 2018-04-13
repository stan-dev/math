#include <stan/math/rev/core/chainablestack.hpp>
#include <iostream>

extern int sassy();

int main() {
  std::cout << "Two's: " << &stan::math::ChainableStack::memalloc_ << std::endl;
  return sassy();
}

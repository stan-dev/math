#include <stan/math/rev/core/allocator.hpp>
#include <stan/math/rev/core/chainablestack.hpp>
#include <vector>
#include <iostream>

using stan::math::ChainableStack;

int main() {
  std::vector<double, autodiff_allocator<double>> ds;
  for (int i = 0; i < 6400; i++) {
    ds.push_back(double(i));
  }

  std::cout << ChainableStack::memalloc_.bytes_allocated() << std::endl;

  for (auto d : ds) {
    d += 1;
  }
}

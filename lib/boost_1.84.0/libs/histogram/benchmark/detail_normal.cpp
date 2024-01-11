// Copyright 2015-2019 Hans Dembinski
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <benchmark/benchmark.h>
#include <boost/histogram/detail/normal.hpp>
#include <boost/math/distributions/normal.hpp>
#include "../test/throw_exception.hpp"
#include "generator.hpp"

namespace bm = boost::math;
using namespace boost::histogram::detail;

#include <cassert>
struct assert_check {
  assert_check() {
    assert(false); // don't run with asserts enabled
  }
} _;

static void math_cdf(benchmark::State& state) {
  generator<uniform> gen;
  bm::normal norm;
  for (auto _ : state) benchmark::DoNotOptimize(cdf(norm, gen()));
  state.SetItemsProcessed(state.iterations());
}

static void our_cdf(benchmark::State& state) {
  generator<uniform> gen;
  for (auto _ : state) benchmark::DoNotOptimize(normal_cdf(gen()));
  state.SetItemsProcessed(state.iterations());
}

static void math_ppf(benchmark::State& state) {
  generator<uniform> gen;
  bm::normal norm;
  for (auto _ : state) benchmark::DoNotOptimize(quantile(norm, gen()));
  state.SetItemsProcessed(state.iterations());
}

static void our_ppf(benchmark::State& state) {
  generator<uniform> gen;
  for (auto _ : state) benchmark::DoNotOptimize(normal_ppf(gen()));
  state.SetItemsProcessed(state.iterations());
}

BENCHMARK(math_cdf);
BENCHMARK(our_cdf);
BENCHMARK(math_ppf);
BENCHMARK(our_ppf);

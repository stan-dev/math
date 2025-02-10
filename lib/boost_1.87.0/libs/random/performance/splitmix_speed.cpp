/*
 * Copyright Matt Borland 2022.
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 */

#include <cstdint>
#include <random>
#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <benchmark/benchmark.h>
#include "../include/boost/random/splitmix64.hpp"

template <typename T>
void stl_mt19937_64(benchmark::State& state)
{
    std::random_device rd;
    std::mt19937_64 gen{rd()};
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
    }
}

template <typename T>
void boost_mt19937_64(benchmark::State& state)
{
    boost::random::random_device rd;
    boost::random::mt19937_64 gen{rd()};
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
        benchmark::DoNotOptimize(gen());
    }
}

template <typename T>
void boost_splitmix64(benchmark::State& state)
{
    boost::random::splitmix64 rd;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(rd());
        benchmark::DoNotOptimize(rd());
        benchmark::DoNotOptimize(rd());
        benchmark::DoNotOptimize(rd());
        benchmark::DoNotOptimize(rd());
        benchmark::DoNotOptimize(rd());
        benchmark::DoNotOptimize(rd());
        benchmark::DoNotOptimize(rd());
        benchmark::DoNotOptimize(rd());
        benchmark::DoNotOptimize(rd());
    }
}

BENCHMARK_TEMPLATE(stl_mt19937_64, std::uint64_t);
BENCHMARK_TEMPLATE(boost_mt19937_64, std::uint64_t);
BENCHMARK_TEMPLATE(boost_splitmix64, std::uint64_t);

BENCHMARK_MAIN();

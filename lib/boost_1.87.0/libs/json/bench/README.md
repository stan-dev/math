# Boost.JSON Benchmarks

To run the benchmarks first run clone.sh to
fetch the third party repositories. Then run
the bench program with no arguments for a
list of command line options.

```
Usage: bench [options...] <file>...

Options:  -t:[p][s]            Test parsing, serialization or both
                                 (default both)
          -i:[b][d][r][c][n]   Test the specified implementations
                                 (b: Boost.JSON, pool storage)
                                 (d: Boost.JSON, default storage)
                                 (u: Boost.JSON, null parser)
                                 (s: Boost.JSON, convenient functions)
                                 (o: Boost.JSON, stream operators)
                                 (r: RapidJSON, memory storage)
                                 (c: RapidJSON, CRT storage)
                                 (n: nlohmann/json)
                                 (default all)
          -n:<number>          Number of trials (default 6)
          -b:<branch>          Branch label for boost implementations
          -m:(i|p|n)           Number parsing mode
                                 (i: imprecise)
                                 (p: precise)
                                 (n: none)
                                 (default imprecise)
          -f                   Include file IO into consideration when testing parsers
```

When building with b2, it is possible to create several different copies of the
bench program for different build properties (toolset, build variant, etc.).
Rather than figuring out where those programs are located from b2 output you
can simply run them directly from the build system. For this, use the
`bench//run` target. For example:

```sh
b2 variant=release toolset=gcc,clang bench//run
```

You can set bench options using the flag `bench.option`. By default the
target uses all files from the `bench/data` directory as benchmark data. If you
want to use just a few of those files or use different files altogether, use
the flag `bench.file`. Finally, when running benchmarks it is often necessary
to use a wrapper program that sets up the process for better result stability.
This can be achieved with the flag `bench.launcher`. In this case also consider
telling b2 to not run several jobs at a time with `-j1`, as they would
interfere with each other. Combining everything together we get:

```sh
b2 -j1 variant=release toolset=gcc,clang bench//run bench.option=-t:p bench.option=-n:10 bench.option=-i:b bench.launcher="nice -n -20 taskset -c 1" bench.file=mydata.json
```

Which is the equivalent of

```sh
nice -n -20 taskset -c 1 path/to/gcc/version/of/bench -t:p -n:10 -i:b mydata.json
nice -n -20 taskset -c 1 path/to/clang/version/of/bench -t:p -n:10 -i:b mydata.json
```

Alternatively, you can copy the version of bench that you want into the bench
directory by requesting the target `bench//bench-local` and then run it
manually:

```sh
b2 variant=release toolset=gcc bench//bench-local
nice -n -20 taskset -c 1 bench/bench-local -t:p -n:10 -i:b mydata.json
```

The benchmarked files were sourced from the
[simdjson](https://github.com/simdjson/simdjson) repository.

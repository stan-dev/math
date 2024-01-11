# Contributing to Boost.Histogram

## Star the project

If you like Boost.Histogram, please star the project on Github! We want Boost.Histogram to be the best histogram library out there. If you give it a star, it becomes more visible and will gain more users. More users mean more user feedback to make the library even better.

## Support

Feel free to ask questions on https://gitter.im/boostorg/histogram.

## Reporting issues

We value your feedback about issues you encounter. The more information you provide the easier it is for developers to resolve the problem.

Issues should be reported to the [issue tracker](
https://github.com/boostorg/histogram/issues?state=open).

Issues can also be used to submit feature requests.

Don't be shy: if you are friendly, we are friendly!

### Reporting failures

Please report any tests failures to the issue tracker along with the test output and information on your system:

* platform (Linux, Windows, OSX, ...)
* compiler and version

## Submitting Pull Requests

Fork the main repository. Base your changes on the `develop` branch. Make a new branch for any feature or bug-fix in your fork. Start developing. You can start a pull request when you feel the change is ready for review.

Please rebase your branch to the original `develop` branch before submitting (which may have diverged from your fork in the meantime).

## Running tests

### With cmake

You can build and test Boost Histogram with `cmake`, but this is **experimental**. It is very convenient, however, as it does not require setting up the Boost superproject, you just need a checkout of Boost Histogram, but when you need to recompile many times, it is slower than using `b2` (see next section). The following steps compile the tests with `cmake`
```sh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
cmake --build . --target check -j3
```
Increase `-j3` to the number of jobs that you want to run in parallel.

### With b2

For general advice on how to set up the Boost project for development, see the Getting Started section on
https://github.com/boostorg/wiki/wiki.

Boost comes with a build system called `b2`, which is the most efficient way to run the develop-test cycle. It takes a few extra steps and some reading to set up, but the payoff is worth it in the long run. After following the general advice from the Boost wiki, you should be all set up to run the tests from the Boost Histogram project folder with the command
```sh
../../b2 cxxstd=latest warnings-as-errors=on test
```
You can also test the examples by executing
```sh
../../b2 cxxstd=latest examples
```
To make the tests complete faster, you can use the option `-j4` (or another number) to build in parallel.

### With b2filt

The output of b2 is verbose and there is no highlighting, so sometimes it is difficult to spot a problem. You can give `b2filt` are try, which you can install with `pip install b2filt`. Just replace `../../b2` with `b2filt`, for example:
```sh
b2filt cxxstd=latest warnings-as-errors=on test
```
Note: b2filt is not an official part of Boost and has never been tested extensively on Windows. Please leave an [issue on GitHub](https://github.com/HDembinski/b2filt) if you have trouble with b2filt.

### Check test coverage

Boost.Histogram maintains 100% line coverage. Coverage is automatically checked by CI. To generate a report locally, you need to build the code with gcc-8 or gcc-12 (these version are known to work, others might) and coverage instrumentation enabled, do `b2 toolset=gcc-8 cxxstd=latest coverage=on test`.

To generate the coverage report, run `tools/cov.py` from the project root directory of Boost.Histogram. This script has not been tested on Windows, so probably will not work on Windows. If you have several versions of gcc installed, you may need to specify the gcov version to use, run for example `GCOV=gcov-8 tools/cov.py`.

This will generate a new folder `coverage-report` with a HTML report. Open `coverage-report/index.html` in a browser.

Notes: Generating coverage data is very fickle. You need to use gcc-8 or gcc-12 and a matching version of gcov. Other gcc versions were known to be broken or were not supported by lcov at the time of this writing, which is used to process the raw coverage data. Generating coverage data with clang and llvm-cov is not supported by lcov at the ime of this writing.

## Generate documentation

To build the documentation, you need to install a few extra things, see [Quickbook documentation](https://www.boost.org/doc/libs/1_74_0/doc/html/quickbook/install.html).

## Coding Style

Follow the [Boost Library Requirements and Guidelines](https://www.boost.org/development/requirements.html) and the established style in Boost.Histogram.

### Code formatting

All names are written with small letters and `_`. Template parameters are capitalized and in camel-case. Do not worry about whitespace, which is automatically formattet with `clang-format`.

To automatically format whitespace in your commits, use `pre-commit`. You can install it with `pip install pre-commit`. Then do `pre-commit install` in the project directory. From now on, `pre-commit` will run clang-format automatically whenever you run `git commit`. If it complains, simply add the changes it made to the files and run `git commit` again.

### Documentation

Doxygen comments should be added for all user-facing functions and methods. Implementation details are not documented (everything in the `boost::histogram::detail` namespace is an implementation detail that can change at any time).

# Contributing to the Stan Math library

Thanks for reading! We love contributions from everyone in the form of good discussion, issues, and pull requests.

This is the short version. For an extended guide see [developer doc on the wiki](https://github.com/stan-dev/math/wiki/Developer-Doc#contributing).

### Pull Request Checklist

1. Fill out a bug or feature request [issue](https://github.com/stan-dev/math/issues).
2. Make a fork and branch, for features start the branch with `feature/your-feature`, for bugs use `bug-fix/your-fix`, and refactors `refactor/your-refactor`.
3. Once your branch passes cpplint and tests relevant to your PR, submit a [pull request](https://github.com/stan-dev/math/pulls) following the checks in the pull request template.

## How to Become a Contributor

### Bug Reports

We reserve [issues](https://github.com/stan-dev/math/issues) for bugs and feature requests that are defined well enough for a developer to tackle. If you have general questions about the Math library, please see the [Discussion](#discussion) section.

The issue template gives the minimum information we like to have when a bug is reported.

1. A description of the problem
2. A reproducible example
3. The expected outcome if the bug were fixed.

If there's an error and you can produce any of these pieces, it's much appreciated!

### Feature Requests

We track the development of new features using the same [issue tracker](https://github.com/stan-dev/math/issues). Ideally, feature requests will have:

1. A description of the feature
2. An example
3. The expected outcome if the feature existed.

Open feature requests should be the ones we want to implement in the Math library. We'll try to close vague feature requests that don't have enough information and move the discussion to the forums.

## Pull Requests

All changes to the Math library are handled through [pull requests](https://github.com/stan-dev/math/pulls). Each pull request should correspond to an issue. We follow a [modified GitFlow branching model](https://github.com/stan-dev/stan/wiki/Dev:-Git-Process) for development.

When a contributor creates a pull request for inclusion to the Math library, here are some of the things we expect:

1. the contribution maintains the Math library's open-source [license](https://github.com/stan-dev/math/wiki/Developer-Doc#licensing): 3-clause BSD
2. the code base remains stable after merging the pull request; we expect the `develop` branch to always be in a good state
3. the changes are maintainable. In code review, we look at the design of the proposed code. We also expect documentation. It should look like idiomatic C++.
4. the changes are tested. For bugs, we expect at least one test that fails before the patch and is fixed after the patch. For new features, we expect at least one test that shows expected behavior and one test that shows the behavior when there's an error.
5. the changes adhere to the Math library's [C++ standards](https://github.com/stan-dev/stan/wiki/Code-Quality). Consistency really helps.

Individual tests can be run with `./runTests.py ./test/path/to/your/tests`. If a path to a folder is given to `./runTests.py` all tests in that folder and subfolders will be run. `make cpplint` will run cpplint to check your changes follow the style guide. We use `clang-format-5.0` to automatically clean formatting. If you have `clang-format-5.0` installed you can run `make clang-format` to run `clang-format` locally, otherwise the Jenkins instance will run `clang-format` and push a commit to your branch with fixes.

TODO: idk where to put the new autodiff testing guide

https://github.com/stan-dev/math/wiki/Automatic-Differentiation-Testing-Framework

After your pull request passes our continuous integration tests, a member of the Stan development team will be assigned to review your code. We expect all the above before a pull request is merged. We are an open-source project and once code makes it into the repository, it's on the community to maintain.

It is the responsibility of the contributor submitting the pull request that the code meets these requirements. We're open-source. Once the code gets into the code base, the community of developers take ownership of it.

### Code Reviews

For a detailed guide See the [Code Review Guidelines](https://github.com/stan-dev/math/wiki/Developer-Doc#code-review-guidelines) on the Math wiki.

For reviewers, please see the [How to Review a Pull Request](https://github.com/stan-dev/math/wiki/Developer-Doc#how-to-review-a-pull-request) in the developer doc.

It's hard to know for a particular PR how much review needs to happen. The below gives some general guidelines, but each reviewer may use their own discretion. Overall, smaller pull requests with clear documentation and tests are the easiest for both the reviewer and reviewee. Medium and large PRs are sometimes able to be broken into multiple smaller PRs. Doing this will almost always lead to an easier review and merge to the math library. All bug fixes and user facing methods and classes must have associated unit tests.


#### Things That Are Easy To Have Accepted

The examples below are generally easy to review pull requests with straightforward testing.

- Bug fixes that have a clear reproducible example
- Small pull requests (< 200 lines of code excluding testing) that add a simple new method or feature
- Documentation
- Anything marked on the issue tracker as `good first issue`

Some easier PRs to include came from below

1. [Rewrite Test Dependencies Script in Python](https://github.com/stan-dev/math/pull/1329)
2. [Adding GP exponential kernel](https://github.com/stan-dev/math/pull/1131/files)
3. [Adding const ref and ref returns to the `to_var` and `to_fvar` methods](https://github.com/stan-dev/math/pull/1323)
4. [Replace hand coded math with C++11 functions](https://github.com/stan-dev/math/pull/1318)
5. [Adding boolean check statements](https://github.com/stan-dev/math/pull/1141)
6. [Support for linear transformed unconstrained parameters](https://github.com/stan-dev/math/pull/1048)

### Things That Are A Little More Difficult

For these you don't need a design document, but your pull request should have a pretty thorough explanation. Something like a small blog post is nice to give the reviewer a bit of background. You may need to give performance tests using the [`perf-math`](https://github.com/seantalts/perf-math) repository.

- Large meaningful refactors
- Big new features

Some example PRs include:

1. [Refactoring type traits](https://github.com/stan-dev/math/pull/1341)
2. [Adding Variadic Templates](https://github.com/stan-dev/math/pull/978)
3. [Allow ts of ODEs to be `var`s](https://github.com/stan-dev/math/pull/857)
4. [Multivariate Container Support For `log_mix`](https://github.com/stan-dev/math/pull/751)
5. [Adding Eigen Plugin Methods for Vars](https://github.com/stan-dev/math/pull/1283)

### Things That Are Going To Be Large Reviews

Think very big PRs or changes to how Stan math fundamentally works. These are a lot easier to review with a very thorough pull request summary and clear tests. Your reviewer may ask you to make a design document which needs to be reviewed before it can be accepted. Depending on the complexity, your reviewer may also request that your large PR be broken up into several smaller PRs.

For some of these we have special cases which I'll list more in depth below

- Changes to how Stan allocates memory
    - Changing how Stan manages memory on the autodiff stack will require performance tests for all compilers and operating systems we support.
- Adding a new dependency
    - See the checklist for [adding a dependency to the math library](https://github.com/stan-dev/math/wiki/Checklist-for-adding-a-dependency-to-the-Math-library) that will need to be filled out.
- Expanding the types of math or models that Stan can work with.
- Changes that heavily effect reproducibility of Stan programs.
- If you need a new folder to hold your code.

These are generally large features which allow Stan to do things it could not previously do, but are also the most complex and difficult to review (for both the reviewer and reviewee). It's heavily encouraged to post on Stan discourse, come to a meeting, and have thorough documentation and discussion before proceeding with a feature heavy pull request.

1. [Async for OpenCL](https://github.com/stan-dev/math/pull/1162)
2. [MPI Backend](https://github.com/stan-dev/math/pull/854)
3. [ODE Solvers](https://github.com/stan-dev/math/pull/768)
4. [Changing how the lgamma function is calculated](https://github.com/stan-dev/math/pull/1255)


## Discussion

For general questions, please ask on the forums with the ["Developers" tag](http://discourse.mc-stan.org/c/stan-dev).


## General Tips

Below are some general tips and good things to know when reading Stan math's code. See the [wiki for new devs](https://github.com/stan-dev/math/wiki/Introduction-to-Stan-Math-for-New-Developers) for an in depth overview for new developers.

### Folder Structure

- `opencl`: The OpenCL backend, everything related to interacting with the OpenCL instance and devices is managed here.
- `memory`: Code for Stan's stack allocator.
- `prim`: Functions here are templated generally such that if a specialization does not exist in rev or fwd the function will default to this implementation in prim. Anything in this folder should be unaware of the `fvar` and `var` types.
- `rev`: Specializations of the functions in prim for reverse mode autodiff.
- `fwd`: Specializations of the functions in prim for forward mode autodiff.
- `mix`: Specializations that operate on mixes of `fvar` and `var` types.

See the [OpenCL Wiki](https://github.com/stan-dev/math/wiki/OpenCL-GPU-Routines) for information on contributing to Stan's GPU backend.

Anything added to `rev` and `fwd` should have a general function available in `prim`. Within each of `prim`, `rev` and `fwd` there are folders for `arr`, `mat` and `scal`. `scal` is for functions that work on scalar types such as `double` or `var`, `arr` holds methods and classes for working with standard library vectors while `mat` holds methods that work on `Eigen` matrices and vectors.


Within each of `scal`, `arr`, and `mat` there are folder for `err`, `fun`, `functor`, `prob` and `meta` which contain

- `err`: Error checking methods.
- `fun`: Functions which operator on types for each folder.
- `functor`: Functors which operate on functions and types for each folder.
- `meta`: Type traits relevant to template metaprogramming.
- `prob`: Methods for calculating the log probability given a distribution.

As it comes to include order, the general rules are

1. `scal` cannot import from `arr` or `mat`.
2. `arr` can include files from `scal`.
3. `mat` can include files from `arr` and `scal`.
4. `prim` cannot include anything from `fwd` or `rev`.
5. `rev` and `fwd` can include from `prim`, but `rev` and `prim` cannot include from one another.
6. `mix` can include from `rev`, `fwd`, and `prim`.

While it makes compile times longer, it's generally best to include the top level file that imports all of the folder. For instance if a method needs to use template metaprogramming methods in any of the `meta` folders for prim, you must include `stan/math/prim/meta.hpp` which imports the meta files in the correct order.

### Adding new methods and using dependencies

It's preferred if methods used within Stan math that are not exposed to the user come from, in order

1. The standard library
2. Boost
3. Eigen and TBB
4. Rolling your own function in the `internal` namespace

Functions in the `internal` namespace should only be used in a single file, else they should be in the `stan math` namespace with proper documentation and tests. It's helpful and good for functions in the `internal` namespace to have at least a one liner comment giving the intent of the function.


Some other random tips

1. Make something `const` if it can be const.
2. Use the necessary error checking available in the `err` folder for your methods inputs and outputs. If your method can't throw an error, make it `noexcept`.

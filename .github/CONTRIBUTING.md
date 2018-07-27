# Contributing to the Stan Math library

Thanks for reading! We love contributions from everyone in the form of good discussion, issues, and pull requests.

This is the short version. There's more information on the [wiki](https://github.com/stan-dev/math/wiki/Developer-Doc#contributing).

## Issues

We reserve [issues](https://github.com/stan-dev/math/issues) for bugs and feature requests that are defined well enough for a developer to tackle. If you have general questions about the Math library, please see the [Discussion](#discussion) section.

### Bug Reports

Ideally, bug reports will include these pieces of information:

1. a description of the problem
2. a reproducible example
3. the expected outcome if the bug were fixed.

If there's an error and you can produce any of these pieces, it's much appreciated!


### Feature Requests

We track the development of new features using the same [issue tracker](https://github.com/stan-dev/math/issues). Ideally, feature requests will have:

1. a description of the feature
2. an example
3. the expected outcome if the feature existed.

Open feature requests should be the ones we want to implement in the Math library. We'll try to close vague feature requests that don't have enough information and move the discussion to the forums.

## Pull Requests

All changes to the Math library are handled through [pull requests](https://github.com/stan-dev/math/pulls). Each pull request should correspond to an issue. We follow a [modified GitFlow branching model](https://github.com/stan-dev/stan/wiki/Dev:-Git-Process) for development.

When a contributor creates a pull request for inclusion to the Math library, here are some of the things we expect:

1. the contribution maintains the Math library's open-source [license](https://github.com/stan-dev/math/wiki/Developer-Doc#licensing): 3-clause BSD
2. the code base remains stable after merging the pull request; we expect the `develop` branch to always be in a good state
3. the changes are maintainable. In code review, we look at the design of the proposed code. We also expect documentation. It should look like idiomatic C++.
4. the changes are tested. For bugs, we expect at least one test that fails before the patch and is fixed after the patch. For new features, we expect at least one test that shows expected behavior and one test that shows the behavior when there's an error.
5. the changes adhere to the Math library's [C++ standards](https://github.com/stan-dev/stan/wiki/Code-Quality). Consistency really helps.

Pull requests are code reviewed after they pass our continuous integration tests. We expect all the above before a pull request is merged. We are an open-source project and once code makes it into the repository, it's on the community to maintain.

It is the responsibility of the contributor submitting the pull request that the code meets these requirements. We're open-source. Once the code gets into the code base, the community of developers take ownership of it.

### Code Reviews
See the [Code Review Guidelines](https://github.com/stan-dev/math/wiki/Developer-Doc#code-review-guidelines) on the Math wiki.


## Discussion

For general questions, please ask on the forums with the ["Developers" tag](http://discourse.mc-stan.org/c/stan-dev).

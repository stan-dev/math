Thanks for submitting a pull request! Please remove this text when submitting.

Start by filling in the Summary, Tests, and Side Effects sections of this pull request and then work through the handy checklist at the bottom. If anything significant is missing, the pull request may be closed until it's ready. The full guidebook on how pull requests are reviewed is here: [Code Review Guidelines](https://github.com/stan-dev/math/wiki/Developer-Doc#code-review-guidelines).

## Summary

Describe the contents of the pull request. The issue describes the problem that needs to be fixed, whereas the text here should help the reviewer figure out to what extent this pull request solves the issue.

Describe implementation details and any relevant information that would help the code reviewer understand this pull request. If there is anything special you need to draw someone's attention to, do it here. Link to any related pull requests or Discourse threads.

## Tests

Describe the new tests with the pull request.

For bug fixes there should be a new test that would fail if the patch weren't in place (so that the bug could be caught if it comes up again).

For new features there should be at least one test showing the expected behavior and one test that demonstrates error handling. Be aware the reviewer will very likely ask for more tests than this, but it's a start.

## Side Effects

Are there any side effects that we should be aware of?

## Release notes

Replace this text with a short note on what will change if this pull request is merged in which case this will be included in the release notes.

## Checklist

- [ ] Math issue #(issue number)

- [ ] Copyright holder: (fill in copyright holder information)

    The copyright holder is typically you or your assignee, such as a university or company. By submitting this pull request, the copyright holder is agreeing to the license the submitted work under the following licenses:
      - Code: BSD 3-clause (https://opensource.org/licenses/BSD-3-Clause)
      - Documentation: CC-BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

- [ ] the basic tests are passing

    - unit tests pass (to run, use: `./runTests.py test/unit`)
    - header checks pass, (`make test-headers`)
    - dependencies checks pass, (`make test-math-dependencies`)
    - docs build, (`make doxygen`)
    - code passes the built in [C++ standards](https://github.com/stan-dev/stan/wiki/Code-Quality) checks (`make cpplint`)

- [ ] the code is written in idiomatic C++ and changes are documented in the doxygen

- [ ] the new changes are tested

Thanks for submitting a pull request! Please remove this text when submitting. We expect everything in the checklist to be completed before the pull request is merged. If any is missing, the pull request may be closed until it's ready. Please help guide the reviewer through your code (see [Code Review Guidelines](https://github.com/stan-dev/math/wiki/Developer-Doc#code-review-guidelines)).

## Checklist

- [ ] Math issue #(issue number)

- [ ] Copyright holder: (fill in copyright holder information)

    The copyright holder is typically you or your assignee, such as a university or company. By submitting this pull request, the copyright holder is agreeing to the license the submitted work under the following licenses:
      - Code: BSD 3-clause (https://opensource.org/licenses/BSD-3-Clause)
      - Documentation: CC-BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

- [ ] the code base is stable

    - all unit tests pass
    - continuous integration passes

- [ ] the changes are maintainable

    - the code design is idiomatic C++ or there's a good reason it's not
    - please include appropriate documentation

- [ ] the changes are tested

    - there's at least one new test
    - for bugs: at least one test that fails before the patch
    - for features: at least one test that shows expected behavior
    - for features: at least one test that shows error handling

- [ ] the changes adhere to the Math library's [C++ standards](https://github.com/stan-dev/stan/wiki/Code-Quality)

    - CppLint passes: `make cpplint`


Additional:
- Link to discourse thread: ()


## Summary
Describe the contents of the pull request. The issue describes what should happen. This pull request should show how that happens.

Describe implementation details and any relevant information that would help the code reviewer understand this pull request.



## Tests
Describe (in words) the new tests with the pull request. Explain what you're testing. There should be at least one test that shows how the new code works. There should also be a test to show what happens when there's a failure.



## Side Effects
Are there any side effects that we should be aware of?


## Additional Notes
Leave additional notes for the code reviewer here.

##
# Stan Math Library
# -----------------
#
# To customize your build, set make variables in either:
#    ~/.config/stan/make.local
#    make/local
# Variables in make/local is loaded after ~/.config/stan/make.local


## 'help' is the default make target.
help:

-include $(HOME)/.config/stan/make.local  # user-defined variables
-include make/local                       # user-defined variables

include make/compiler_flags               # CXX, CXXFLAGS, LDFLAGS set by the end of this file
include make/dependencies                 # rules for generating dependencies
include make/libraries
include make/tests
include make/cpplint

.PHONY: help
help:
	@echo '--------------------------------------------------------------------------------'
	@echo 'Note: testing of Math is typically done with the `runTests.py` python script.'
	@echo '  See https://github.com/stan-dev/math/wiki/Developer-Doc#building-and-running-tests'
	@echo '  for more detail on testing.'
	@echo  ''
	@echo 'Stan Math makefile:'
	@$(MAKE) print-compiler-flags
	@echo 'Tests:'
	@echo ''
	@echo '  Unit tests are built through make by specifying the executable as the target'
	@echo '  to make. For a test in test/*_test.cpp, the executable is test/*$(EXE).'
	@echo ''
	@echo '  Header tests'
	@echo '  - test-headers  : tests all source headers to ensure they are compilable and'
	@echo '                    include enough header files.'
	@echo ''
	@echo '  To run a single header test, add "-test" to the end of the file name.'
	@echo '  Example: make stan/math/constants.hpp-test'
	@echo ''
	@echo '  - test-math-dependencies : walks through all the header files and indicates'
	@echo '      when the math dependencies are violated. Dependencies should follow:'
	@echo '      * rev -> prim'
	@echo '      * fwd -> prim'
	@echo '      * mix -> {rev, fwd, prim}'
	@echo '      * within {prim, rev, fwd, mix}: mat -> arr -> scal'
	@echo ''
	@echo '  Cpplint'
	@echo '  - cpplint       : runs cpplint.py on source files. requires python 2.7.'
	@echo '                    cpplint is called using the CPPLINT variable:'
	@echo '                      CPPLINT = $(CPPLINT)'
	@echo '                    To set the version of python 2, set the PYTHON2 variable:'
	@echo '                      PYTHON2 = $(PYTHON2)'
	@echo ''
	@echo 'Documentation:'
	@echo '  Doxygen'
	@echo '  - doxygen       : runs doxygen on source files. requires doxygen.'
	@echo ''
	@echo 'Clean:'
	@echo '  - clean         : Basic clean. Leaves doc and compiled libraries intact.'
	@echo '  - clean-deps    : Removes dependency files for tests. If tests stop building,'
	@echo '                    run this target.'
	@echo '  - clean-libraries : Removes binaries built for libraries including CVODES.'
	@echo '  - clean-all     : Cleans up all of Stan.'
	@echo ''
	@echo '--------------------------------------------------------------------------------'


.PHONY: doxygen
doxygen:
	mkdir -p doc/api
	doxygen doxygen/doxygen.cfg

##
# Clean up.
##
.PHONY: clean clean-doxygen clean-deps clean-all
clean:
	@echo '  removing test executables'
	$(shell find test -type f -name "*_test$(EXE)" -exec rm {} +)
	$(shell find test -type f -name "*_test.d" -exec rm {} +)
	$(shell find test -type f -name "*_test.d.*" -exec rm {} +)
	$(shell find test -type f -name "*_test.xml" -exec rm {} +)
	$(shell find test -type f -name "*.o" -exec rm {} +)
	$(RM) $(wildcard $(GTEST)/src/gtest-all.o)
	@echo '  removing generated test files'
	$(RM) $(wildcard test/prob/generate_tests$(EXE))
	$(shell find test/prob -name '*_generated_*_test.cpp' -type f -exec rm {} +)

clean-doxygen:
	@echo '  removing doxygen'
	$(RM) -r doc/api

clean-deps:
	@echo '  removing dependency files'
	$(shell find stan test lib -type f -name '*.d' -exec rm {} +)
	$(shell find stan test lib -type f -name '*.d.*' -exec rm {} +)
	$(shell find stan  -type f -name '*.dSYM' -exec rm {} +)

clean-all: clean clean-doxygen clean-deps clean-libraries

##
# Debug target that allows you to print a variable
##
print-%  : ; @echo $* = $($*)

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
include make/clang-tidy

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
	@echo '      * only include {prim, rev, fwd, mix}/meta.hpp from the meta subfolders'
	@echo ''
	@echo '  Cpplint'
	@echo '  - cpplint       : runs cpplint.py on source files. requires python 2.7.'
	@echo '                    cpplint is called using the CPPLINT variable:'
	@echo '                      CPPLINT = $(CPPLINT)'
	@echo '                    To set the version of python 2, set the PYTHON2 variable:'
	@echo '                      PYTHON2 = $(PYTHON2)'
	@echo ''
	@echo ' Clang Tidy'
	@echo ' - clang-tidy     : runs the clang-tidy makefile over the test suite.'
	@echo '                    Options:'
	@echo '                     files: (Optional) regex for file names to include in the check'
	@echo '                      Default runs all the tests in unit'
	@echo '                     tidy_checks: (Optional) A set of checks'
	@echo '                      Default runs a hand picked selection of tests'
	@echo ''
	@echo '     Example: This runs clang-tidy over all the multiply tests in prim'
	@echo ''
	@echo '     make clang-tidy files=*prim*multiply*'
	@echo ''
	@echo ' - clang-tidy-fix : same as above but runs with the -fix flag.'
	@echo '                    For automated fixes, outputs a yaml named'
	@echo '                    .clang-fixes.yml'
	@echo ''
	@echo 'Documentation:'
	@echo '  Doxygen'
	@echo '  - doxygen       : runs doxygen on source files. requires doxygen.'
	@echo ''
	@echo 'Clean:'
	@echo '  - clean         : Basic clean. Leaves doc and compiled libraries intact.'
	@echo '  - clean-deps    : Removes dependency files for tests. If tests stop building,'
	@echo '                    run this target.'
	@echo '  - clean-libraries : Removes binaries built for libraries including CVODES and the TBB.'
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
	@$(RM) $(call findfiles,test,*_test$(EXE))
	@$(RM) $(call findfiles,test,*_test.d)
	@$(RM) $(call findfiles,test,*_test.d.*)
	@$(RM) $(call findfiles,test,*_test.xml)
	@$(RM) $(call findfiles,test,*.o)
	@$(RM) $(wildcard $(GTEST)/src/gtest-all.o)
	@echo '  removing generated test files'
	@$(RM) $(wildcard test/prob/generate_tests$(EXE))
	@$(RM) $(call findfiles,test/prob,*_generated_*_test.cpp)

clean-doxygen:
	@echo '  removing doxygen'
	$(RM) -r doc/api

clean-deps:
	@echo '  removing dependency files'
	@$(RM) $(call findfiles,stan,*.d)
	@$(RM) $(call findfiles,test,*.d)
	@$(RM) $(call findfiles,lib,*.d)
	@$(RM) $(call findfiles,stan,*.d.*)
	@$(RM) $(call findfiles,test,*.d.*)
	@$(RM) $(call findfiles,lib,*.d.*)
	@$(RM) $(call findfiles,stan,*.dSYM)

clean-all: clean clean-doxygen clean-deps clean-libraries

.PHONY: test-math-dependencies
test-math-dependencies:
	@./runChecks.py
##
# Debug target that allows you to print a variable
##
print-%  : ; @echo $* = $($*)

##
# Stan Math Library
# -----------------
#
# To customize your build, set make variables in either:
#    ~/.config/stan/make.local
#    make/local
# Variables in make/local is loaded after ~/.config/stan/make.local
##

# 'help' is the default make target.
help:

-include $(HOME)/.config/stan/make.local  # user-defined variables
-include make/local                       # user-defined variables

include make/defaults
include make/libraries
include make/tests
include make/cpplint



.PHONY: help
help:
	@echo '--------------------------------------------------------------------------------'
	@echo 'Stan Math Library'
	@echo '  https://github.com/stan-dev/math'
	@echo ''
	@echo '  Current configuration:'
	@echo '  - CXX                            ' $(CXX)
	@echo '  - O (optimization level):        ' $(O)
	@echo '  - CXXFLAGS:                      ' $(CXXFLAGS)
	@echo '  - AR (archiver):                 ' $(AR)
	@echo '  - OS                             ' $(OS)
	@echo ''
	@echo '  Library configuration:'
	@echo '  - EIGEN                       ' $(EIGEN)
	@echo '  - BOOST                       ' $(BOOST)
	@echo '  - IDAS                        ' $(IDAS)
	@echo '  - CVODES                      ' $(CVODES)
	@echo '  - GTEST                       ' $(GTEST)
	@echo ''
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

## doxygen
.PHONY: doxygen
doxygen:
	mkdir -p doc/api
	doxygen doxygen/doxygen.cfg

##
# Clean up.
##
.PHONY: clean clean-doxygen clean-deps clean-all
clean:
	@echo 'removing test executables'
	$(shell find test -type f -name "*_test$(EXE)" -exec rm {} +)
	$(shell find test -type f -name "*_test.d" -exec rm {} +)
	$(shell find test -type f -name "*_test.d.*" -exec rm {} +)
	$(shell find test -type f -name "*_test.xml" -exec rm {} +)
	$(shell find test -type f -name "*.o" -exec rm {} +)
	$(shell find test -type f -name "lib*.so" -exec rm {} +)

clean-doxygen:
	$(RM) -r doc/api

clean-deps:
	@echo 'removing dependency files'
	$(shell find . -type f -name '*.d' -exec rm {} +)
	$(shell find . -type f -name '*.d.*' -exec rm {} +)
	$(RM) $(shell find stan -type f -name '*.dSYM') $(shell find stan -type f -name '*.d.*')

clean-all: clean clean-doxygen clean-deps clean-libraries
	@echo 'removing generated test files'
	$(shell find test/prob -name '*_generated_*_test.cpp' -type f -exec rm {} +)
	$(RM) $(wildcard test/prob/generate_tests$(EXE))

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





# ##
# # Get information about the compiler used.
# # - CC_TYPE: {g++, clang++, mingw32-g++, other}
# # - CC_MAJOR: major version of CC
# # - CC_MINOR: minor version of CC
# ##
# -include make/detect_cc

# # OS_TYPE is set automatically by this script
# ##
# # These includes should update the following variables
# # based on the OS:
# #   - CFLAGS
# #   - GTEST_CXXFLAGS
# #   - EXE
# ##
# -include make/detect_os

# include make/tests    # tests
# include make/cpplint  # cpplint

# ##
# # Dependencies
# ##
# ifneq (,$(filter-out test-headers generate-tests clean% %-test %.d,$(MAKECMDGOALS)))
#   -include $(addsuffix .d,$(subst $(EXE),,$(MAKECMDGOALS)))
# endif


.PHONY: help
help:
	@echo 'MAKECMDGOALS: ' $(MAKECMDGOALS)
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
# 	@echo '  - OS_TYPE (Operating System): ' $(OS_TYPE)
# 	@echo '  - CC (Compiler):              ' $(CC)
# 	@echo '  - CC_TYPE                     ' $(CC_TYPE)
# 	@echo '  - Compiler version:           ' $(CC_MAJOR).$(CC_MINOR)
# 	@echo '  - O (Optimization Level):     ' $(O)
# 	@echo '  Library configuration:'
# 	@echo '  - EIGEN                       ' $(EIGEN)
# 	@echo '  - BOOST                       ' $(BOOST)
# 	@echo '  - CVODES                      ' $(CVODES)
# 	@echo '  - GTEST                       ' $(GTEST)
# 	@echo ''
# 	@echo 'Tests:'
# 	@echo ''
# 	@echo '  Unit tests are built through make by specifying the executable as the target'
# 	@echo '  to make. For a test in test/*_test.cpp, the executable is test/*$(EXE).'
# 	@echo ''
# 	@echo '  Header tests'
# 	@echo '  - test-headers  : tests all source headers to ensure they are compilable and'
# 	@echo '                    include enough header files.'
# 	@echo ''
# 	@echo '  To run a single header test, add "-test" to the end of the file name.'
# 	@echo '  Example: make stan/math/constants.hpp-test'
# 	@echo ''
# 	@echo '  - test-math-dependencies : walks through all the header files and indicates'
# 	@echo '      when the math dependencies are violated. Dependencies should follow:'
# 	@echo '      * rev -> prim'
# 	@echo '      * fwd -> prim'
# 	@echo '      * mix -> {rev, fwd, prim}'
# 	@echo '      * within {prim, rev, fwd, mix}: mat -> arr -> scal'
# 	@echo ''
# 	@echo '  Cpplint'
# 	@echo '  - cpplint       : runs cpplint.py on source files. requires python 2.7.'
# 	@echo '                    cpplint is called using the CPPLINT variable:'
# 	@echo '                      CPPLINT = $(CPPLINT)'
# 	@echo '                    To set the version of python 2, set the PYTHON2 variable:'
# 	@echo '                      PYTHON2 = $(PYTHON2)'
# 	@echo ''
# 	@echo 'Documentation:'
# 	@echo '  Doxygen'
# 	@echo '  - doxygen       : runs doxygen on source files. requires doxygen.'
# 	@echo ''
# 	@echo 'Clean:'
# 	@echo '  - clean         : Basic clean. Leaves doc and compiled libraries intact.'
# 	@echo '  - clean-deps    : Removes dependency files for tests. If tests stop building,'
# 	@echo '                    run this target.'
# 	@echo '  - clean-libraries : Removes binaries built for libraries including CVODES.'
# 	@echo '  - clean-all     : Cleans up all of Stan.'
# 	@echo ''
	@echo '--------------------------------------------------------------------------------'

# ## doxygen
# .PHONY: doxygen
# doxygen:
# 	mkdir -p doc/api
# 	doxygen doxygen/doxygen.cfg

# ##
# # Clean up.
# ##
# .PHONY: clean clean-doxygen clean-deps clean-all
# clean:
# 	@echo '  removing test executables'
# 	$(shell find test -type f -name "*_test$(EXE)" -exec rm {} +)
# 	$(shell find test -type f -name "*_test.d" -exec rm {} +)
# 	$(shell find test -type f -name "*_test.d.*" -exec rm {} +)
# 	$(shell find test -type f -name "*_test.xml" -exec rm {} +)
# 	$(shell find test -type f -name "*.o" -exec rm {} +)
# 	$(shell find test -type f -name "lib*.so" -exec rm {} +)

# clean-doxygen:
# 	$(RM) -r doc/api

# clean-deps:
# 	@echo '  removing dependency files'
# 	$(shell find . -type f -name '*.d' -exec rm {} +)
# 	$(shell find . -type f -name '*.d.*' -exec rm {} +)
# 	$(RM) $(shell find stan -type f -name '*.dSYM') $(shell find stan -type f -name '*.d.*')

# clean-all: clean clean-doxygen clean-deps clean-libraries
# 	@echo '  removing generated test files'
# 	$(shell find test/prob -name '*_generated_*_test.cpp' -type f -exec rm {} +)
# 	$(RM) $(wildcard test/prob/generate_tests$(EXE))

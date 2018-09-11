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
	@echo 'Stan Math makefile:'
	@echo '  Current configuration:'
	@echo '  - OS (Operating System):      ' $(OS)
	@echo '  - CXX (Compiler):             ' $(CXX)
	@echo '  - CXX_TYPE                    ' $(CXX_TYPE)
	@echo '  - Compiler version:           ' $(CXX_MAJOR).$(CXX_MINOR)
	@echo '  - O (Optimization Level):     ' $(O)
	@echo '  Library configuration:'
	@echo '  - EIGEN                       ' $(EIGEN)
	@echo '  - BOOST                       ' $(BOOST)
	@echo '  - SUNDIALS                    ' $(SUNDIALS)
	@echo '  - GTEST                       ' $(GTEST)
	@echo '  - STAN_OPENCL                 ' $(STAN_OPENCL)
	@echo '  - STAN_MPI                    ' $(STAN_MPI)
	@echo '  Compiler flags (each can be overriden separately):'
	@echo '  - CXXFLAGS_LANG               ' $(CXXFLAGS_LANG)
	@echo '  - CXXFLAGS_WARNINGS           ' $(CXXFLAGS_WARNINGS)
	@echo '  - CXXFLAGS_BOOST              ' $(CXXFLAGS_BOOST)
	@echo '  - CXXFLAGS_EIGEN              ' $(CXXFLAGS_EIGEN)
	@echo '  - CXXFLAGS_OS                 ' $(CXXFLAGS_OS)
	@echo '  - CXXFLAGS_GTEST              ' $(CXXFLAGS_GTEST)
	@echo '  - CXXFLAGS_OPENCL             ' $(CXXFLAGS_OPENCL)
	@echo '  - CXXFLAGS_MPI                ' $(CXXFLAGS_MPI)
	@echo '  - CFLAGS_SUNDIALS             ' $(CFLAGS_SUNDIALS)
	@echo '  LDLIBS:'
	@echo '  - LDLIBS_LANG                 ' $(LDLIBS_LANG)
	@echo '  - LDLIBS_WARNINGS             ' $(LDLIBS_WARNINGS)
	@echo '  - LDLIBS_BOOST                ' $(LDLIBS_BOOST)
	@echo '  - LDLIBS_EIGEN                ' $(LDLIBS_EIGEN)
	@echo '  - LDLIBS_OS                   ' $(LDLIBS_OS)
	@echo '  - LDLIBS_GTEST                ' $(LDLIBS_GTEST)
	@echo '  - LDLIBS_OPENCL               ' $(LDLIBS_OPENCL)
	@echo '  - LDLIBS_MPI                  ' $(LDLIBS_MPI)
	@echo '  LDFLAGS:'
	@echo '  - LDFLAGS_LANG                ' $(LDFLAGS_LANG)
	@echo '  - LDFLAGS_WARNINGS            ' $(LDFLAGS_WARNINGS)
	@echo '  - LDFLAGS_BOOST               ' $(LDFLAGS_BOOST)
	@echo '  - LDFLAGS_EIGEN               ' $(LDFLAGS_EIGEN)
	@echo '  - LDFLAGS_OS                  ' $(LDFLAGS_OS)
	@echo '  - LDFLAGS_GTEST               ' $(LDFLAGS_GTEST)
	@echo '  - LDFLAGS_OPENCL              ' $(LDFLAGS_OPENCL)
	@echo '  - LDFLAGS_MPI                 ' $(LDFLAGS_MPI)
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

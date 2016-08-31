#include <gtest/gtest.h>

TEST(multiple_translation_units, compile) {
  SUCCEED()
    << "this test compiling indicates that compiling the math library "
    << "with multiple translation units is ok.";
}

"""Exercise the production CMake test helpers with small compile/link fixtures."""

import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent


@unittest.skipUnless(shutil.which("cmake"), "CMake is required")
class CMakeBuildTest(unittest.TestCase):
    def check_build(self, pch_enabled):
        with tempfile.TemporaryDirectory(prefix="stan-cmake-build-") as temporary:
            source = Path(temporary) / "source"
            build = Path(temporary) / "build"
            commands = Path(temporary) / "compile.jsonl"

            def write(path, text):
                destination = source / path
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_text(text)

            def run(*args):
                result = subprocess.run(args, text=True, stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT)
                self.assertEqual(result.returncode, 0, result.stdout)

            write("launcher.py", """import json, subprocess, sys
with open(sys.argv[1], 'a') as log:
    log.write(json.dumps(sys.argv[2:]) + '\\n')
sys.exit(subprocess.call(sys.argv[2:]))
""")
            write("CMakeLists.txt", f"""
cmake_minimum_required(VERSION 3.16.3)
project(stan_cmake_fixture LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(STAN_TEST_PCH {'ON' if pch_enabled else 'OFF'})
set(STAN_TEST_OPTIMIZATION_FLAGS -O1)
add_library(stanmath INTERFACE)
target_include_directories(stanmath INTERFACE "${{CMAKE_CURRENT_SOURCE_DIR}}")
add_library(gtest_main STATIC main.cpp)
enable_testing()
include("{(ROOT / 'test/CMakeLists.txt').as_posix()}")
""")
            write("cmake/pch.cpp", "// Empty executable/PCH source.\n")
            write("gtest/gtest.h", """#pragma once

#include <vector>

using test_callback = bool (*)();
inline std::vector<test_callback>& test_callbacks() {
  static std::vector<test_callback> callbacks;
  return callbacks;
}
struct test_registration {
  explicit test_registration(test_callback callback) {
    test_callbacks().push_back(callback);
  }
};
""")
            write("main.cpp", """#include <gtest/gtest.h>
int main() {
  for (auto callback : test_callbacks()) {
    if (!callback()) return 1;
  }
  return 0;
}
""")
            write("stan/math/prim.hpp", """#pragma once

#ifdef STAN_PCH_TEST_SETTING
constexpr int configured_value = STAN_PCH_TEST_SETTING;
#else
constexpr int configured_value = 1;
#endif
""")
            for layer in ("fwd", "rev", "mix"):
                write(f"stan/math/{layer}.hpp", "#include <stan/math/prim.hpp>\n")
            write("stan/math/mix.hpp", "#include <stan/math/prim.hpp>\nusing mix_only_alias = double;\n")
            write("test/unit/math/test_ad.hpp", """#pragma once

#include <stan/math/mix.hpp>
#include <gtest/gtest.h>

using test_alias = double;
""")
            write("unit/math/prim/fun/plain_test.cpp", """#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

using mix_only_alias = int;
using test_alias = int;
bool plain_test() { return configured_value == 1; }
test_registration plain_registration(plain_test);
""")
            write("unit/math/prim/fun/ad_test.cpp", """#include <test/unit/math/test_ad.hpp>
bool ad_test() { return configured_value == 1 && sizeof(test_alias) == sizeof(double); }
test_registration ad_registration(ad_test);
""")
            write("unit/math/prim/fun/configuration_test.cpp", """#define STAN_PCH_TEST_SETTING 2

#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

static_assert(configured_value == 2, "Configuration must precede Stan headers");
bool configuration_test() { return configured_value == 2; }
test_registration configuration_registration(configuration_test);
""")
            write("generated_template.cpp", """#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

bool generated_test() { return configured_value == 1; }
test_registration generated_registration(generated_test);
""")
            write("unit/CMakeLists.txt", """
foreach(layer prim fwd rev mix)
  stan_add_pch(unit_${layer}_pch <gtest/gtest.h> <stan/math/${layer}.hpp>)
endforeach()
stan_add_pch(unit_ad_pch <test/unit/math/test_ad.hpp>)
add_gtest_grouped_test("${CMAKE_CURRENT_SOURCE_DIR}/math" unit_math unit_mix_pch)
set(generated "${PROJECT_BINARY_DIR}/generated/generated_test.cpp")
add_custom_command(OUTPUT "${generated}"
  COMMAND "${CMAKE_COMMAND}" -E make_directory "${PROJECT_BINARY_DIR}/generated"
  COMMAND "${CMAKE_COMMAND}" -E copy "${PROJECT_SOURCE_DIR}/generated_template.cpp" "${generated}"
  DEPENDS "${PROJECT_SOURCE_DIR}/generated_template.cpp")
add_custom_target(generate_fixture DEPENDS "${generated}")
set_source_files_properties("${generated}" PROPERTIES STAN_GENERATOR_TARGET generate_fixture)
add_gtest_grouped_test("${PROJECT_BINARY_DIR}/generated" expressions unit_prim_pch "${generated}")
""")
            launcher = f"{sys.executable};{source / 'launcher.py'};{commands}"
            run("cmake", "-S", str(source), "-B", str(build),
                "-DCMAKE_CXX_COMPILER_LAUNCHER=" + launcher)
            run("cmake", "--build", str(build), "--parallel", "1", "--target",
                "unit_math_prim_fun_all_tests", "expressions_all_tests")
            run(str(build / "unit/unit_math_prim_fun_all_tests"))
            run(str(build / "unit/expressions_all_tests"))

            compiled = [json.loads(line) for line in commands.read_text().splitlines()]
            fixtures = {name: source / f"unit/math/prim/fun/{name}_test.cpp"
                        for name in ("plain", "ad", "configuration")}
            fixtures["generated"] = build / "generated/generated_test.cpp"
            expected_pch = {"plain": "unit_prim_pch", "ad": "unit_ad_pch",
                            "configuration": None, "generated": "unit_prim_pch"}
            for name, path in fixtures.items():
                invocations = [args for args in compiled if str(path) in args]
                self.assertEqual(len(invocations), 1, (name, invocations))
                pch = expected_pch[name] if pch_enabled else None
                pch_arguments = [arg for arg in invocations[0] if "cmake_pch" in arg]
                if pch:
                    self.assertTrue(any(pch in arg for arg in pch_arguments), invocations[0])
                else:
                    self.assertFalse(pch_arguments, invocations[0])

            # Switching executables must not compile any test source again.
            run("cmake", "--build", str(build), "--parallel", "1", "--target",
                "test_unit_math_prim_fun_plain_test", "unit_math_all_tests",
                "test_expressions_generated_test")
            run(str(build / "unit/test_unit_math_prim_fun_plain_test"))
            run(str(build / "unit/unit_math_all_tests"))
            run(str(build / "unit/test_expressions_generated_test"))
            after = [json.loads(line) for line in commands.read_text().splitlines()]
            for name, path in fixtures.items():
                self.assertEqual(sum(str(path) in args for args in after), 1, name)

    def test_pch_profiles_opt_out_generation_and_object_reuse(self):
        self.check_build(True)

    def test_without_pch_preserves_object_reuse(self):
        self.check_build(False)


if __name__ == "__main__":
    unittest.main()

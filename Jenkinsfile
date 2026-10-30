properties([
  disableConcurrentBuilds(abortPrevious: true),
  buildDiscarder(logRotator(numToKeepStr: '20', daysToKeepStr: '30')),
  parameters([
    string(defaultValue: 'develop', name: 'cmdstan_pr', description: 'PR to test CmdStan upstream against e.g. PR-630'),
    string(defaultValue: 'develop', name: 'stan_pr', description: 'PR to test Stan upstream against e.g. PR-630'),
    booleanParam(defaultValue: false, name: 'withRowVector', description: 'Run additional distribution tests on RowVectors (takes 5x as long)'),
    booleanParam(defaultValue: false, name: 'disableJumbo', description: 'Disable Jumbo tests. This takes longer and should only be used for debugging if it is believed that the jumbo tests are causing failures.'),
    booleanParam(defaultValue: false, name: 'optimizeUnitTests', description: 'Use O=3 for unit tests (takex ~3x as long)'),
    booleanParam(defaultValue: false, name: 'runAllDistributions', description: 'Run all distribution tests, even ones which are unchanged compared to develop'),
    booleanParam(defaultValue: false, name: 'run_all', description: 'Pretend all files changes'),
  ])
])

def image = 'stanorg/ci:v1'
def runRemainingStages = false
def CLANG_CXX = 'clang++-7'
def CLANG_CC = 'clang-7'
def GCC = 'g++'
def MPICXX = 'mpicxx.openmpi'

def mainBranch = env.BRANCH_NAME == 'develop' || env.BRANCH_NAME == 'master'
def noOptimize = !(params.optimizeUnitTests || mainBranch)

// Quote complete arguments, including compiler flags containing spaces.
def shellQuote(String value) {
  return "'" + value.replace("'", "'\"'\"'") + "'"
}

// Every configure of `build` starts from a fresh cache, so it uses exactly
// these settings. stanorg/ci:v1 has CMake 3.16, which lacks `cmake --fresh`
// (3.24), so the cache and top-level CMakeFiles are deleted instead, which is
// what --fresh does. Objects whose flags are unchanged are not rebuilt.
def cmakeSettings(List settings) {
  return ['-DSTAN_USE_MAKE_LOCAL=OFF'] + settings
}

def configureCmake(List settings) {
  sh 'cmake -E remove -f build/CMakeCache.txt'
  sh 'cmake -E remove_directory build/CMakeFiles'
  sh 'cmake -S . -B build ' + cmakeSettings(settings).collect { shellQuote(it.toString()) }.join(' ')
}

// runTests.py --fresh does the same cache reset before it configures.
def runnerSettings(List settings) {
  return (['--build-dir', 'build', '--fresh'] + cmakeSettings(settings).collect {
    "--cmake-arg=${it}"
  }).collect { shellQuote(it.toString()) }.join(' ')
}

def runTests(List settings, String args) {
  sh 'cmake -E remove_directory build/test-results'
  sh "python3 runTests.py -j\$PARALLEL ${runnerSettings(settings)} $args"
  junit(testResults: 'build/test-results/*.xml', allowEmptyResults: true)
}

catchError {
  withEnv([
    'GCC=g++',
    'GIT_AUTHOR_NAME=Stan Jenkins',
    'GIT_AUTHOR_EMAIL=mc.stanislaw@gmail.com',
    'GIT_COMMITTER_NAME=Stan Jenkins',
    'GIT_COMMITTER_EMAIL=mc.stanislaw@gmail.com'
  ]) {
    runPod(image: image, cpus: 2) {
      stage('Verify changes') {
        runRemainingStages = params.run_all || filesChanged(
            'stan', 'make', 'lib', 'test', 'cmake', 'CMakeLists.txt', 'runTests.py', 'runChecks.py', 'makefile', 'Jenkinsfile', '.clang-format')
      }

      stage("Clang-format") {
        def dirty = sh returnStatus: true, script: """
          clang-format --version
          git ls-files 'stan/*.hpp' 'test/*.hpp' 'stan/*.cpp' 'test/*.cpp' | xargs -n20 -P\$PARALLEL clang-format -i
          git diff --exit-code
        """
        if (dirty) {
          def branch = env.CHANGE_BRANCH ?: env.BRANCH_NAME
          def repo = env.CHANGE_FORK ?: "stan-dev"
          if (!("/" in repo))
            repo += "/math.git"
          echo "Exiting build because clang-format found changes."
          emailext (
              subject: "[StanJenkins] Autoformattted: Job '${env.JOB_NAME} [${env.BUILD_NUMBER}]'",
              body: """
Job '${env.JOB_NAME} [${env.BUILD_NUMBER}]' has been autoformatted and the
changes committed to your branch, if permissions allowed.  Please pull these
changes before continuing.

See https://github.com/stan-dev/stan/wiki/Coding-Style-and-Idioms for setting
up the autoformatter locally.  (Check console output at ${env.BUILD_URL})
""",
              recipientProviders: [[$class: 'RequesterRecipientProvider']],
              to: env.CHANGE_AUTHOR_EMAIL)
          sh '''
            git add -u stan test
            git commit -m "[Jenkins] auto-formatting by `clang-format --version`"
          '''
          gitPush(gitScm: scmGit(
              userRemoteConfigs: [[credentialsId: "stan-github", name: 'dest', url: "https://github.com/$repo"]],
              branches: [[name: "refs/heads/$branch"]]),
              targetBranch: branch,
              targetRepo: 'dest')
          echo "Those changes are now found on stan-dev/math under $repo branch $branch"
          echo "Please 'git pull' before continuing to develop."
          error "clang-format changes"
        }
      }

      stage('Linting & Doc checks') {
        configureCmake(["-DCMAKE_C_COMPILER=${CLANG_CC}", "-DCMAKE_CXX_COMPILER=${CLANG_CXX}",
                        '-DBUILD_TESTING=OFF', '-DSTAN_BUILD_DOCS=ON'])
        parallel(
          CppLint: { sh "cmake --build build --target cpplint" },
          Dependencies: { sh "cmake --build build --target test-math-dependencies" },
          Documentation: { sh "cmake --build build --target doxygen" },
        )
      }
    }

    if (runRemainingStages) {
      stage('Quick tests') {
        parallel headers: {
          runPod(image: image) {
            stage('Headers check') {
              configureCmake(["-DCMAKE_C_COMPILER=${CLANG_CC}", "-DCMAKE_CXX_COMPILER=${CLANG_CXX}",
                              '-DCMAKE_CXX_FLAGS=-Werror', '-DSTAN_TEST_HEADERS=ON',
                              '-DSTAN_BUILD_PROB_TESTS=OFF'])
              sh 'cmake --build build --parallel $PARALLEL --target test-headers'
            }
          }
        }, unit: {
          if (env.CHANGE_TARGET) {
            runPod(image: image) {
              stage('Run changed unit tests') {
                runTests(['-DSTAN_TEST_OPTIMIZATION=3', '-DCMAKE_CXX_FLAGS=-march=native -mtune=native'],
                         '--changed --debug')
              }
            }
          }
        }
      }

      stage('Full Unit Tests') {
        parallel failFast: true,
          revfwd: {
            runPod(image: image) {
              stage('Rev/Fwd Unit Tests') {
                def settings = ['-DSTAN_SANITIZE_ADDRESS=ON',
                                "-DSTAN_TEST_OPTIMIZATION=${noOptimize ? '0' : '3'}"]
                runTests(settings, "test/unit/math/rev")
                runTests(settings, "test/unit/math/fwd")
              }
            }
          },
          mix: {
            runPod(image: image, cpus: 8, memory: '192Gi') {
              stage('Mix Unit Tests') {
                def settings = ['-DSTAN_SANITIZE_ADDRESS=ON',
                                "-DSTAN_TEST_OPTIMIZATION=${noOptimize ? '1' : '3'}"]
                // A single mix executable exceeds the 2 GiB x86-64 small code
                // model limit under ASan, so link each subfolder separately.
                runTests(settings, "test/unit/math/mix/*_test.cpp")
                runTests(settings, "test/unit/math/mix/constraint")
                runTests(settings, "test/unit/math/mix/core")
                runTests(settings, "test/unit/math/mix/fun")
                runTests(settings, "test/unit/math/mix/functor")
                runTests(settings, "test/unit/math/mix/meta")
                runTests(settings, "test/unit/math/mix/prob")
              }
            }
          },
          prim: {
            runPod(image: image) {
              stage('Prim Unit Tests') {
                def settings = ['-DSTAN_SANITIZE_ADDRESS=ON',
                                "-DSTAN_TEST_OPTIMIZATION=${noOptimize ? '0' : '3'}"]
                runTests(settings, "test/unit/*_test.cpp")
                runTests(settings, "test/unit/math/*_test.cpp")
                runTests(settings, "test/unit/math/prim")
                runTests(settings, "test/unit/math/memory")
              }
            }
          },
          laplace: {
            runPod(image: image, memory: '32Gi') {
              stage('Laplace Unit Tests') {
                def settings = ['-DCMAKE_CXX_FLAGS=-march=native -mtune=native', '-DSTAN_TEST_OPTIMIZATION=3']

                runTests(settings, "test/unit/math/laplace/*_test.cpp")

                if (mainBranch) {
                  // only run ASAN on a selected test to avoid extremely long CI times
                  runTests(settings + ['-DSTAN_SANITIZE_ADDRESS=ON'],
                           "test/unit/math/laplace/laplace_marginal_lpdf_moto_test.cpp")
                }
              }
            }
          },
          gpu: {
            runPod(image: image, memory: '32Gi', gpus: 1) {
              stage('OpenCL GPU tests') {
                def settings = ["-DCMAKE_C_COMPILER=${CLANG_CC}", "-DCMAKE_CXX_COMPILER=${CLANG_CXX}",
                                '-DCMAKE_CXX_FLAGS=-Werror', '-DSTAN_OPENCL=ON',
                                '-DSTAN_OPENCL_PLATFORM_ID=0', '-DSTAN_OPENCL_DEVICE_ID=0',
                                '-DOpenCL_LIBRARY=/usr/local/cuda/targets/x86_64-linux/lib/libOpenCL.so',
                                "-DSTAN_TEST_OPTIMIZATION=${noOptimize ? '1' : '3'}"]
                runTests(settings, "test/unit/math/opencl") // TODO(bward): try to enable
                runTests(settings, "test/unit/multiple_translation_units_test.cpp")
              }
            }
          },
          mpi: {
            runPod(image: image) {
              stage('MPI Tests') {
                def settings = ["-DCMAKE_CXX_COMPILER=${MPICXX}", '-DSTAN_MPI=ON', '-DSTAN_TEST_OPTIMIZATION=3']
                runTests(settings, "test/unit/math/prim/functor")
                runTests(settings, "test/unit/math/rev/functor")
              }
            }
          },
          expr: {
            runPod(image: image, memory: '64Gi') {
              stage('Expressions test') {
                def settings = ["-DCMAKE_C_COMPILER=${CLANG_CC}", "-DCMAKE_CXX_COMPILER=${CLANG_CXX}",
                                '-DCMAKE_CXX_FLAGS=-Werror', '-DSTAN_TEST_OPTIMIZATION=0']
                sh "python ./test/code_generator_test.py"
                sh "python ./test/signature_parser_test.py"
                sh "python ./test/statement_types_test.py"
                sh "python ./test/varmat_compatibility_summary_test.py"
                // These scripts configure their own CMake builds; CC, CXX, and
                // CXXFLAGS set those builds' compilers and CMAKE_CXX_FLAGS.
                withEnv(["CC=${CLANG_CC}", "CXX=${CLANG_CXX}", 'CXXFLAGS=-Werror']) {
                  sh "python ./test/varmat_compatibility_test.py"
                  sh "python ./test/expressions/test_expression_testing_framework.py"
                }
                runTests(settings, 'test/expressions')
                runTests(settings + ['-DSTAN_THREADS=ON'], 'test/expressions --only-functions reduce_sum map_rect')
              }
            }
          },
          thread: {
            runPod(image: image, cpus: 8, memory: '128Gi') {
              stage('Threading tests') {
                def settings = ["-DCMAKE_C_COMPILER=${CLANG_CC}", "-DCMAKE_CXX_COMPILER=${CLANG_CXX}",
                                '-DCMAKE_CXX_FLAGS=-Werror', '-DSTAN_THREADS=ON', '-DSTAN_TEST_OPTIMIZATION=3']
                withEnv(['STAN_NUM_THREADS=4', 'PARALLEL=8']) { // runTests currently only runs 1 test at a time post-build
                  if (mainBranch) {
                    runTests(settings, "test/unit")
                  } else {
                    runTests(settings, "test/unit -f thread")
                    runTests(settings, "test/unit -f map_rect")
                    runTests(settings, "test/unit -f reduce_sum")
                  }
                }
              }
            }
          },
          dist: {
          def cores_per_test = 6
          def parallel_tests = 10
          runPod(image: image, cpus: cores_per_test*parallel_tests, memory: '128Gi') {
            stage ('Distribution tests') {
              def settings = ["-DCMAKE_C_COMPILER=${CLANG_CC}", "-DCMAKE_CXX_COMPILER=${CLANG_CXX}",
                              '-DSTAN_TEST_OPTIMIZATION=0', '-DSTAN_PROB_TEST_BATCH_SIZE=100']
              if (params.withRowVector || mainBranch) {
                settings += ['-DSTAN_TEST_ROW_VECTORS=ON', '-DSTAN_PROB_TEST_ALL=ON']
              }
              // getDependencies.py lists each fixture's includes with the
              // configured build's compiler flags.
              sh "python3 runTests.py --configure-only ${runnerSettings(settings)}"

              def cmd = 'python3 test/prob/getDependencies.py --build-dir build'
              if (params.runAllDistributions || mainBranch) {
                cmd += ' --pretend-all'
              }
              sh """
              $cmd > changed-tests
              if [ -s changed-tests ] ; then
                ./runTests.py -j${cores_per_test*parallel_tests} ${runnerSettings(settings)} --make-only `cat changed-tests`
                unset PARALLEL
                parallel --halt now,fail=1 -r -j$parallel_tests ./runTests.py --build-dir build --test-only -j$cores_per_test {} < changed-tests
              fi
              """
            }
          }
        }
      }

      if (env.CHANGE_TARGET) {
        stage('Upstream tests') {
          build(job: "CCM/Stan/stan/$params.stan_pr",
            parameters: [
              booleanParam(name: 'downstream', value: true),
              string(name: 'math_pr', value: env.BRANCH_NAME),
              string(name: 'cmdstan_pr', value: params.cmdstan_pr)
            ])
        }
      }
    }

    if (env.BRANCH_NAME == 'develop') {
      runPod(image: image, cpus: 2) {
        stage('Upload doxygen') {
          def branch = 'gh-pages'
          configureCmake(['-DBUILD_TESTING=OFF', '-DSTAN_BUILD_DOCS=ON'])
          sh 'cmake --build build --target doxygen'
          dir('gh-pages') {
            def scm = scmGit(
              branches: [[name: "refs/heads/$branch"]],
              userRemoteConfigs: [[credentialsId: 'stan-github', url: 'https://github.com/stan-dev/math.git']])
            checkout(scm: scm, changelog: false, poll: false)
            sh '''
              rm -rf *
              mv ../doc/api/html/* .
              git add .
              git commit --amend -m "auto generated docs from Jenkins"
            '''
            withCredentials([gitUsernamePassword(credentialsId: 'stan-github', gitToolName: 'git-tool')]) {
              sh "git push --force origin HEAD:$branch"
            }
          }
        }

        stage('Update upstream') {
          dir('stan') {
            updateSubmodule('stan', 'develop', 'lib/stan_math')
          }
        }
      }
    }
  }
}

emailFailure()

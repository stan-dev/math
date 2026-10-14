properties([
  disableConcurrentBuilds(abortPrevious: true),
  buildDiscarder(logRotator(numToKeepStr: '20', daysToKeepStr: '30')),
  parameters([
    string(defaultValue: 'develop', name: 'cmdstan_pr', description: 'PR to test CmdStan upstream against e.g. PR-630'),
    string(defaultValue: 'develop', name: 'stan_pr', description: 'PR to test Stan upstream against e.g. PR-630'),
    booleanParam(defaultValue: false, name: 'withRowVector', description: 'Run additional distribution tests on RowVectors (takes 5x as long)'),
    booleanParam(defaultValue: false, name: 'optimizeUnitTests', description: 'Use O=3 for unit tests (takes ~3x as long)'),
    booleanParam(defaultValue: false, name: 'runAllDistributions', description: 'Run all distribution tests, even ones which are unchanged compared to develop'),
    booleanParam(defaultValue: false, name: 'run_all', description: 'Pretend all files changed'),
  ])
])

def image = 'stanorg/ci:gpu-cmake3.30.5'
def runRemainingStages = false
def CLANG_CXX = 'clang++'
def CLANG_CC = 'clang'
def MPICXX = 'mpicxx.openmpi'

def mainBranch = env.BRANCH_NAME == 'develop' || env.BRANCH_NAME == 'master'
def noOptimize = !(params.optimizeUnitTests || mainBranch)

// Quote complete arguments, including regexes and compiler flags containing spaces.
def shellQuote(String value) {
  return "'" + value.replace("'", "'\"'\"'") + "'"
}

def runCmakeTests(Map config) {
  def buildDir = config.buildDir
  def jobs = Math.max(1, Math.min(6, (env.PARALLEL ?: '6').toInteger()))
  def flags = [
    '-S', '.', '-B', buildDir, '-DCMAKE_BUILD_TYPE=Release',
    '-DSTAN_TEST_HEADERS=OFF', '-DSTAN_BUILD_PROB_TESTS=OFF',
    "-DSTAN_TEST_FILTER=${config.filter}",
    "-DSTAN_TEST_OPTIMIZATION=${config.optimization ?: '1'}"
  ] + (config.flags ?: [])
  try {
    sh "cmake -E rm -rf ${shellQuote(buildDir + '/test-results')} " +
       "${shellQuote(buildDir + '/ctest-results.xml')} ${shellQuote(buildDir + '/Testing')}"
    sh "cmake " + flags.collect { shellQuote(it.toString()) }.join(' ')
    sh "cmake --build ${shellQuote(buildDir)} --parallel ${jobs}"
    sh "cmake -E make_directory ${shellQuote(buildDir + '/test-results')}"
    withEnv(["GTEST_OUTPUT=xml:${pwd()}/${buildDir}/test-results/"]) {
      sh "ctest --test-dir ${shellQuote(buildDir)} --parallel ${config.testJobs ?: jobs} " +
         "--no-tests=error --output-on-failure " +
         "--output-junit ${shellQuote(pwd() + '/' + buildDir + '/ctest-results.xml')}"
    }
  } finally {
    junit(testResults: "${buildDir}/test-results/*.xml", allowEmptyResults: true)
    archiveArtifacts(artifacts: "${buildDir}/ctest-results.xml,${buildDir}/Testing/Temporary/LastTest.log",
                     allowEmptyArchive: true)
  }
}

// Var-matrix compatibility still invokes Make directly; expression framework
// self-tests exercise compile and runtime failures through the CMake runner.
def runLegacyFrameworkTests(String compiler) {
  writeFile(file: 'make/local', text: "CXX=${compiler}\nO=0\n")
  withEnv(['PARALLEL=6', 'MAKEFLAGS=-j1', 'PATH+TBB=./lib/tbb', "CXX=${compiler}"]) {
    sh 'python3 test/varmat_compatibility_test.py'
    sh 'python3 test/expressions/test_expression_testing_framework.py'
  }
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
        writeFile(file: "make/local", text: "CXX=$CLANG_CXX\nBOOST_PARALLEL_JOBS=\$PARALLEL")
        parallel(
          CppLint: { sh "make cpplint" },
          Dependencies: { sh "make test-math-dependencies" },
          CMakeSelection: {
            sh 'python3 test/cmake_ci_test.py'
            sh 'python3 test/cmake_runner_test.py'
            sh 'python3 test/cmake_build_test.py'
          },
          Documentation: {
            sh 'cmake -S . -B build-ci-docs -DBUILD_TESTING=OFF -DSTAN_BUILD_DOCS=ON'
            sh 'cmake --build build-ci-docs --parallel 6 --target doxygen'
          },
        )
      }
    }

    if (runRemainingStages) {
      stage('Quick tests') {
        parallel headers: {
          runPod(image: image) {
            stage('Headers check') {
              sh "cmake -S . -B build-ci-headers -DCMAKE_BUILD_TYPE=Release " +
                 "-DCMAKE_C_COMPILER=${CLANG_CC} -DCMAKE_CXX_COMPILER=${CLANG_CXX} " +
                 "-DCMAKE_CXX_FLAGS=-Werror -DSTAN_BUILD_PROB_TESTS=OFF " +
                 '-DSTAN_TEST_FILTER=\'^$\' -DSTAN_TEST_HEADERS=ON'
              sh 'cmake --build build-ci-headers --parallel 6 --target test-headers'
            }
          }
        }, unit: {
          if (env.CHANGE_TARGET) {
            runPod(image: image) {
              stage('Run changed unit tests') {
                def filter = sh(returnStdout: true,
                  script: 'python3 test/cmake_ci.py changed-unit-filter').trim()
                if (filter) {
                  runCmakeTests(buildDir: 'build-ci-changed', filter: filter, optimization: '3',
                    flags: ['-DCMAKE_CXX_FLAGS=-march=native -mtune=native'])
                } else {
                  echo 'No priority unit tests selected.'
                }
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
                runCmakeTests(buildDir: 'build-ci-revfwd', filter: '^test/unit/math/(rev|fwd)/',
                  optimization: noOptimize ? '0' : '3',
                  flags: ['-DCMAKE_CXX_FLAGS=-fsanitize=address', '-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address',
                          '-DCMAKE_SHARED_LINKER_FLAGS=-fsanitize=address'])
              }
            }
          },
          mix: {
            runPod(image: image, cpus: 6, memory: '192Gi') {
              stage('Mix Unit Tests') {
                runCmakeTests(buildDir: 'build-ci-mix', filter: '^test/unit/math/mix/',
                  optimization: noOptimize ? '1' : '3',
                  flags: ['-DCMAKE_CXX_FLAGS=-fsanitize=address', '-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address',
                          '-DCMAKE_SHARED_LINKER_FLAGS=-fsanitize=address'])
              }
            }
          },
          prim: {
            runPod(image: image) {
              stage('Prim Unit Tests') {
                runCmakeTests(buildDir: 'build-ci-prim',
                  filter: '^test/unit/([^/]+_test[.]cpp$|math/([^/]+_test[.]cpp$|prim/|memory/))',
                  optimization: noOptimize ? '0' : '3',
                  flags: ['-DCMAKE_CXX_FLAGS=-fsanitize=address', '-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address',
                          '-DCMAKE_SHARED_LINKER_FLAGS=-fsanitize=address'])
              }
            }
          },
          laplace: {
            runPod(image: image, memory: '32Gi') {
              stage('Laplace Unit Tests') {
                runCmakeTests(buildDir: 'build-ci-laplace', filter: '^test/unit/math/laplace/',
                  optimization: '3', flags: ['-DCMAKE_CXX_FLAGS=-march=native -mtune=native'])
                if (mainBranch) {
                  runCmakeTests(buildDir: 'build-ci-laplace-asan', optimization: '3',
                    filter: '^test/unit/math/laplace/laplace_marginal_lpdf_moto_test[.]cpp$',
                    flags: ['-DCMAKE_CXX_FLAGS=-march=native -mtune=native -fsanitize=address',
                            '-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address'])
                }
              }
            }
          },
          gpu: {
            runPod(image: image, memory: '32Gi', gpus: 1) {
              stage('OpenCL GPU tests') {
                runCmakeTests(buildDir: 'build-ci-opencl',
                  filter: '^test/unit/(math/opencl/|multiple_translation_units_test[.]cpp$)',
                  optimization: noOptimize ? '1' : '3',
                  flags: ["-DCMAKE_C_COMPILER=${CLANG_CC}", "-DCMAKE_CXX_COMPILER=${CLANG_CXX}",
                          '-DCMAKE_CXX_FLAGS=-Werror', '-DSTAN_OPENCL=ON',
                          '-DSTAN_OPENCL_PLATFORM_ID=0', '-DSTAN_OPENCL_DEVICE_ID=0',
                          '-DOpenCL_LIBRARY=/usr/local/cuda/targets/x86_64-linux/lib/libOpenCL.so'])
              }
            }
          },
          mpi: {
            runPod(image: image) {
              stage('MPI Tests') {
                runCmakeTests(buildDir: 'build-ci-mpi',
                  filter: '^test/unit/math/(prim|rev)/functor/', testJobs: 1,
                  flags: ["-DCMAKE_CXX_COMPILER=${MPICXX}", '-DSTAN_MPI=ON'])
              }
            }
          },
          expr: {
            runPod(image: image, memory: '64Gi') {
              stage('Expressions test') {
                sh 'python3 test/code_generator_test.py'
                sh 'python3 test/signature_parser_test.py'
                sh 'python3 test/statement_types_test.py'
                sh 'python3 test/varmat_compatibility_summary_test.py'
                runLegacyFrameworkTests(CLANG_CXX)
                sh 'python3 test/cmake_ci.py generate-expressions build-ci-expressions'
                runCmakeTests(buildDir: 'build-ci-expressions', filter: '^test/expressions/', optimization: '0',
                  flags: ["-DCMAKE_C_COMPILER=${CLANG_CC}", "-DCMAKE_CXX_COMPILER=${CLANG_CXX}",
                          '-DCMAKE_CXX_FLAGS=-Werror', '-DSTAN_BUILD_EXPRESSION_TESTS=ON'])
                sh 'python3 test/cmake_ci.py generate-expressions build-ci-expressions-threads --functions reduce_sum map_rect'
                withEnv(['STAN_NUM_THREADS=4']) {
                  runCmakeTests(buildDir: 'build-ci-expressions-threads', filter: '^test/expressions/',
                    optimization: '0', testJobs: 1,
                    flags: ["-DCMAKE_C_COMPILER=${CLANG_CC}", "-DCMAKE_CXX_COMPILER=${CLANG_CXX}",
                            '-DCMAKE_CXX_FLAGS=-Werror', '-DSTAN_BUILD_EXPRESSION_TESTS=ON', '-DSTAN_THREADS=ON'])
                }
              }
            }
          },
          thread: {
            runPod(image: image, cpus: 6, memory: '128Gi') {
              stage('Threading tests') {
                withEnv(['STAN_NUM_THREADS=4']) {
                  runCmakeTests(buildDir: 'build-ci-threads', testJobs: 1,
                    filter: mainBranch ? '^test/unit/' : '^test/unit/.*(thread|map_rect|reduce_sum)',
                    flags: ["-DCMAKE_C_COMPILER=${CLANG_CC}", "-DCMAKE_CXX_COMPILER=${CLANG_CXX}",
                            '-DCMAKE_CXX_FLAGS=-Werror', '-DSTAN_THREADS=ON'])
                }
              }
            }
          },
          dist: {
            runPod(image: image, cpus: 6, memory: '128Gi') {
              stage('Distribution tests') {
                def command = 'python3 test/cmake_ci.py changed-distribution-filter'
                if (params.runAllDistributions || params.run_all || mainBranch) {
                  command += ' --all'
                }
                def filter = sh(returnStdout: true, script: command).trim()
                if (filter) {
                  def flags = ["-DCMAKE_C_COMPILER=${CLANG_CC}", "-DCMAKE_CXX_COMPILER=${CLANG_CXX}",
                               '-DSTAN_BUILD_PROB_TESTS=ON', '-DSTAN_PROB_TEST_BATCH_SIZE=100']
                  if (params.withRowVector || mainBranch) {
                    flags += ['-DSTAN_TEST_ROW_VECTORS=ON', '-DSTAN_PROB_TEST_ALL=ON']
                  }
                  runCmakeTests(buildDir: 'build-ci-distributions', filter: filter,
                    optimization: '0', flags: flags)
                } else {
                  echo 'No distribution tests affected by this change.'
                }
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
          sh 'cmake -S . -B build-ci-docs -DBUILD_TESTING=OFF -DSTAN_BUILD_DOCS=ON'
          sh 'cmake --build build-ci-docs --parallel 6 --target doxygen'
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

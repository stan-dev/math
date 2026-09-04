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
def GCC = 'g++'
def MPICXX = 'mpicxx.openmpi'

def mainBranch = env.BRANCH_NAME == 'develop' || env.BRANCH_NAME == 'master'
def noOptimize = !(params.optimizeUnitTests || mainBranch)
def jumboFlags = params.disableJumbo ? '' : ' --jumbo --debug'

def runTests(String local, String args) {
  catchError(buildResult: 'UNSTABLE', stageResult: 'UNSTABLE') {
    writeFile(file: "make/local", text: local)
    sh "cat make/local"
    sh "make print-compiler-flags"
    sh "python3 runTests.py -j\$PARALLEL $args"
  }
  junit 'test/**/*.xml'
  sh "find test -name *.xml -delete"
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
            'stan', 'make', 'lib', 'test', 'runTests.py', 'runChecks.py', 'makefile', 'Jenkinsfile', '.clang-format')
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
            git add -u src
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
          Documentation: { sh "make doxygen" },
        )
        /* TODO: recordIsuses? */
      }
    }

    if (runRemainingStages) {
      stage('Quick tests') {
        parallel headers: {
          runPod(image: image) {
            stage('Headers check') {
              writeFile(file: 'make/local', text: "CXX=$CLANG_CXX -Werror")
              sh 'make -j$PARALLEL test-headers'
            }
          }
        }, unit: {
          if (env.CHANGE_TARGET) {
            runPod(image: image) {
              stage('Run changed unit tests') {
                writeFile(file: 'make/local', text: "O=3\nCXXFLAGS+=-march=native -mtune=native")
                sh './runTests.py -j$PARALLEL --changed --debug'
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
                def local = 'CXXFLAGS+= -fsanitize=address\n'
                if (noOptimize)
                  local += 'O=0\n'
                runTests(local, "test/unit/math/rev")
                runTests(local, "test/unit/math/fwd")
              }
            }
          },
          mix: {
            runPod(image: image, cpus: 8, memory: '192Gi') {
              stage('Mix Unit Tests') {
                def local = 'CXXFLAGS+= -fsanitize=address\n'
                if (noOptimize)
                  local += 'O=1\n'
                runTests(local, "test/unit/math/mix" + jumboFlags)
              }
            }
          },
          prim: {
            runPod(image: image) {
              stage('Prim Unit Tests') {
                def local = 'CXXFLAGS+= -fsanitize=address\n'
                if (noOptimize)
                  local += 'O=0\n'
                runTests(local, "test/unit/*_test.cpp")
                runTests(local, "test/unit/math/*_test.cpp")
                runTests(local, "test/unit/math/prim" + jumboFlags)
                runTests(local, "test/unit/math/memory")
              }
            }
          },
          laplace: {
            runPod(image: image, memory: '32Gi') {
              stage('Laplace Unit Tests') {
                def local = 'CXXFLAGS+= -march=native -mtune=native\nO=3\n'
                if (!noOptimize)
                  local += 'CXXFLAGS+= -fsanitize=address'
                runTests(local, "test/unit/math/laplace/*_test.cpp")
              }
            }
          },
          gpu: {
            runPod(image: image, memory: '32Gi', gpus: 1) {
              stage('OpenCL GPU tests') {
                def local = """CXX=$CLANG_CXX -Werror
STAN_OPENCL=true
OPENCL_PLATFORM_ID_GPU=0
OPENCL_DEVICE_ID=0
LDFLAGS_OPENCL=-L/usr/local/cuda/targets/x86_64-linux/lib
"""
                if (noOptimize)
                  local += 'O=1\n'
                runTests(local, "test/unit/math/opencl") // TODO(bward): try to enable
                runTests(local, "test/unit/multiple_translation_units_test.cpp")
              }
            }
          }
      }

      stage('Always-run tests') {
        parallel failFast: true,
          mpi: {
            runPod(image: image) {
              stage('Laplace Unit Tests') {
                def local = "CXX=$MPICXX\nCXX_TYPE=gcc\nSTAN_MPI=true\n"
                runTests(local, "test/unit/math/prim/functor")
                runTests(local, "test/unit/math/rev/functor")
              }
            }
          },
          expr: {
            runPod(image: image, memory: '64Gi') {
              stage('Expressions test') {
                def local = "CXX=$CLANG_CXX -Werror\nO=0\n"
                writeFile(file: 'make/local', text: local)
                sh "python ./test/code_generator_test.py"
                sh "python ./test/signature_parser_test.py"
                sh "python ./test/statement_types_test.py"
                sh "python ./test/varmat_compatibility_summary_test.py"
                sh "python ./test/varmat_compatibility_test.py"
                withEnv(['PATH+TBB=./lib/tbb']) {
                  sh "python ./test/expressions/test_expression_testing_framework.py"
                  runTests(local, 'test/expressions')
                  sh "make clean-all"
                  runTests(local + "STAN_THREADS=true\n", 'test/expressions --only-functions reduce_sum map_rect')
                }
              }
            }
          },
          thread: {
            runPod(image: image, memory: '32Gi') {
              stage('Threading tests') {
                def local = "CXX=$CLANG_CXX -Werror\nSTAN_THREADS=true\n"
                if (mainBranch) {
                  runTests(local, "test/unit")
                } else {
                  runTests(local, "test/unit -f thread")
                  runTests(local, "test/unit -f map_rect")
                  runTests(local, "test/unit -f reduce_sum")
                }
              }
            }
          }
      }

      def cores_per_test = 6
      def parallel_tests = 10
      runPod(image: image, cpus: cores_per_test*parallel_tests, memory: '128Gi') {
        stage ('Distribution tests') {
          def local = "CXX=$CLANG_CXX\nO=0\nN_TESTS=100\n"
          if (params.withRowVector || mainBranch) {
            local += "CXXFLAGS+= -DSTAN_TEST_ROW_VECTORS -DSTAN_PROB_TEST_ALL\n"
          }
          writeFile(file: 'make/local', text: local)

          def cmd = 'python3 test/prob/getDependencies.py'
          if (params.runAllDistributions || mainBranch) {
            cmd += ' --pretend-all'
          }
          sh """
            $cmd > changed-tests
            if [ -s changed-tests ] ; then
              ./runTests.py -j${cores_per_test*parallel_tests} --make-only `cat changed-tests`
              unset PARALLEL
              parallel --halt now,fail=1 -r -j$parallel_tests ./runTests.py --test-only -j$cores_per_test {} < changed-tests
            fi
          """
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
          sh 'make doxygen'
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

#!/usr/bin/env groovy
@Library('StanUtils')
import org.stan.Utils

def runTests(String testPath, boolean jumbo = false) {
    try {
        if (jumbo && !params.disableJumbo) {
            sh "python3 runTests.py -j${env.PARALLEL} ${testPath} --jumbo --debug"
        } else {
            sh "python3 runTests.py -j${env.PARALLEL} ${testPath}"
        }
    }
        finally { junit 'test/**/*.xml' }
}

def skipRemainingStages = false
def changedDistributionTests = []
def utils = new org.stan.Utils()

def isBranch(String b) { env.BRANCH_NAME == b }

String alsoNotify() {
    if (isBranch('master') || isBranch('develop')) {
        "stan-buildbot@googlegroups.com"
    } else ""
}
Boolean isPR() { env.CHANGE_URL != null }
String fork() { env.CHANGE_FORK ?: "stan-dev" }
String branchName() { isPR() ? env.CHANGE_BRANCH :env.BRANCH_NAME }
String cmdstan_pr() { params.cmdstan_pr ?: ( env.CHANGE_TARGET == "master" ? "downstream_hotfix" : "downstream_tests" ) }
String stan_pr() { params.stan_pr ?: ( env.CHANGE_TARGET == "master" ? "downstream_hotfix" : "downstream_tests" ) }

pipeline {
    agent none
    parameters {
        string(defaultValue: '', name: 'cmdstan_pr', description: 'PR to test CmdStan upstream against e.g. PR-630')
        string(defaultValue: '', name: 'stan_pr', description: 'PR to test Stan upstream against e.g. PR-630')
        booleanParam(defaultValue: false, name: 'withRowVector', description: 'Run additional distribution tests on RowVectors (takes 5x as long)')
        booleanParam(defaultValue: false, name: 'disableJumbo', description: 'Disable Jumbo tests. This takes longer and should only be used for debugging if it is believed that the jumbo tests are causing failures.')
        booleanParam(defaultValue: false, name: 'optimizeUnitTests', description: 'Use O=3 for unit tests (takex ~3x as long)')
        booleanParam(defaultValue: false, name: 'runAllDistributions', description: 'Run all distribution tests, even ones which are unchanged compared to develop')
    }
    options {
        skipDefaultCheckout()
        preserveStashes(buildCount: 7)
        parallelsAlwaysFailFast()
    }
    environment {
        STAN_NUM_THREADS = 4
        CLANG_CXX = 'clang++-7'
        GCC = 'g++'
        MPICXX = 'mpicxx.openmpi'
        N_TESTS = 100
        OPENCL_DEVICE_ID = 0
        OPENCL_DEVICE_ID_CPU = 0
        OPENCL_DEVICE_ID_GPU = 0
        OPENCL_PLATFORM_ID = 1
        OPENCL_PLATFORM_ID_CPU = 0
        OPENCL_PLATFORM_ID_GPU = 0
        PARALLEL = 4
        GIT_AUTHOR_NAME = 'Stan Jenkins'
        GIT_AUTHOR_EMAIL = 'mc.stanislaw@gmail.com'
        GIT_COMMITTER_NAME = 'Stan Jenkins'
        GIT_COMMITTER_EMAIL = 'mc.stanislaw@gmail.com'
    }
    stages {

        stage('Kill previous builds') {
            when {
                not { branch 'develop' }
                not { branch 'master' }
            }
            steps {
                script {
                    utils.killOldBuilds()
                }
            }
        }

        stage("Clang-format") {
            agent {
                docker {
                    image 'stanorg/ci:gpu-cpp17'
                    label 'linux'
                }
            }
            steps {
                retry(3) { checkout scm }
                withCredentials([usernamePassword(credentialsId: 'a630aebc-6861-4e69-b497-fd7f496ec46b',
                    usernameVariable: 'GIT_USERNAME', passwordVariable: 'GIT_PASSWORD')]) {
                    sh """#!/bin/bash
                        set -x
                        git checkout -b ${branchName()}
                        clang-format --version
                        find stan test -name '*.hpp' -o -name '*.cpp' | xargs -n20 -P${PARALLEL} clang-format -i
                        if [[ `git diff` != "" ]]; then
                            git add stan test
                            git commit -m "[Jenkins] auto-formatting by `clang-format --version`"
                            git push https://${GIT_USERNAME}:${GIT_PASSWORD}@github.com/${fork()}/math.git ${branchName()}
                            echo "Exiting build because clang-format found changes."
                            echo "Those changes are now found on stan-dev/math under branch ${branchName()}"
                            echo "Please 'git pull' before continuing to develop."
                            exit 1
                        fi"""
                }
            }
            post {
                always { deleteDir() }
                failure {
                    script {
                        emailext (
                            subject: "[StanJenkins] Autoformattted: Job '${env.JOB_NAME} [${env.BUILD_NUMBER}]'",
                            body: "Job '${env.JOB_NAME} [${env.BUILD_NUMBER}]' " +
                                "has been autoformatted and the changes committed " +
                                "to your branch, if permissions allowed." +
                                "Please pull these changes before continuing." +
                                "\n\n" +
                                "See https://github.com/stan-dev/stan/wiki/Coding-Style-and-Idioms" +
                                " for setting up the autoformatter locally.\n"+
                            "(Check console output at ${env.BUILD_URL})",
                            recipientProviders: [[$class: 'RequesterRecipientProvider']],
                            to: "${env.CHANGE_AUTHOR_EMAIL}"
                        )
                    }
                }
            }
        }

        stage('Linting & Doc checks') {
            agent {
                docker {
                    image 'stanorg/ci:gpu-cpp17'
                    label 'linux'
                }
            }
            steps {
                script {
                    retry(3) { checkout scm }
                    sh "git clean -xffd"
                    stash 'MathSetup'
                    sh "echo CXX=${CLANG_CXX} > make/local"
                    sh "echo BOOST_PARALLEL_JOBS=${PARALLEL} >> make/local"
                    parallel(
                        CppLint: { sh "make cpplint" },
                        Dependencies: { sh """#!/bin/bash
                            set -o pipefail
                            make test-math-dependencies 2>&1 | tee dependencies.log""" } ,
                        Documentation: { sh "make doxygen" },
                    )
                }
            }
            post {
                always {
                    recordIssues enabledForFailure: true, tools:
                        [cppLint(),
                         groovyScript(parserId: 'mathDependencies', pattern: '**/dependencies.log')]
                    deleteDir()

                }
            }
        }

        stage('Verify changes') {
            agent {
                docker {
                    image 'stanorg/ci:gpu-cpp17'
                    label 'linux'
                }
            }
            steps {
                script {

                    retry(3) { checkout scm }
                    sh 'git clean -xffd'

                    def paths = ['stan', 'make', 'lib', 'test', 'runTests.py', 'runChecks.py', 'makefile', 'Jenkinsfile', '.clang-format'].join(" ")
                    skipRemainingStages = utils.verifyChanges(paths)

                }
            }
            post { always { retry(3) { deleteDir() } } }
        }

        stage('Quick tests') {
            when {
                expression {
                    !skipRemainingStages
                }
            }
            failFast true
            parallel {
                stage('Headers check') {
                    when {
                        expression {
                            !skipRemainingStages
                        }
                    }
                    agent {
                        docker {
                            image 'stanorg/ci:gpu-cpp17'
                            label 'linux'
                        }
                    }

                    steps {
                        unstash 'MathSetup'
                        sh "echo CXX=${CLANG_CXX} -Werror > make/local"
                        sh "make -j${PARALLEL} test-headers"
                    }
                    post { always { deleteDir() } }
                }
                stage('Run changed unit tests') {
                    agent {
                        docker {
                            image 'stanorg/ci:gpu-cpp17'
                            label 'linux'
                            args '--cap-add SYS_PTRACE'
                        }
                    }
                    when {
                        allOf {
                            expression {
                                env.BRANCH_NAME ==~ /PR-\d+/
                            }
                            expression {
                                !skipRemainingStages
                            }
                        }
                    }
                    steps {
                        retry(3) { checkout scm }

                        sh "echo CXXFLAGS += -fsanitize=address >> make/local"
                        sh "./runTests.py -j${PARALLEL} --changed --debug"

                    }
                    post { always { retry(3) { deleteDir() } } }
                }
            }
        }

        stage('Full Unit Tests') {
            when {
                expression {
                    !skipRemainingStages
                }
            }
            failFast true
            parallel {
                stage('Rev/Fwd Unit Tests') {
                    agent {
                        docker {
                            image 'stanorg/ci:gpu-cpp17'
                            label 'linux'
                            args '--cap-add SYS_PTRACE'
                        }
                    }
                    when {
                        expression {
                            !skipRemainingStages
                        }
                    }
                    steps {
                        unstash 'MathSetup'
                        sh "echo CXXFLAGS += -fsanitize=address >> make/local"

                        script {
                            if (!(params.optimizeUnitTests || isBranch('develop') || isBranch('master'))) {
                                sh "echo O=0 >> make/local"
                            }

                            runTests("test/unit/math/rev")
                            runTests("test/unit/math/fwd")
                        }
                    }
                    post { always { retry(3) { deleteDir() } } }
                }
                stage('Mix Unit Tests') {
                    agent {
                        docker {
                            image 'stanorg/ci:gpu-cpp17'
                            label 'linux'
                            args '--cap-add SYS_PTRACE'
                        }
                    }
                    when {
                        expression {
                            !skipRemainingStages
                        }
                    }
                    steps {
                        unstash 'MathSetup'
                        sh "echo CXXFLAGS += -fsanitize=address >> make/local"
                        script {
                            if (!(params.optimizeUnitTests || isBranch('develop') || isBranch('master'))) {
                                sh "echo O=1 >> make/local"
                            }
                            runTests("test/unit/math/mix", true)
                        }
                    }
                    post { always { retry(3) { deleteDir() } } }
                }
                stage('Prim Unit Tests') {
                    agent {
                        docker {
                            image 'stanorg/ci:gpu-cpp17'
                            label 'linux'
                            args '--cap-add SYS_PTRACE'
                        }
                    }
                    when {
                        expression {
                            !skipRemainingStages
                        }
                    }
                    steps {
                        unstash 'MathSetup'
                        sh "echo CXXFLAGS += -fsanitize=address >> make/local"
                        script {
                            if (!(params.optimizeUnitTests || isBranch('develop') || isBranch('master'))) {
                                sh "echo O=0 >> make/local"
                            }
                            runTests("test/unit/*_test.cpp", false)
                            runTests("test/unit/math/*_test.cpp", false)
                            runTests("test/unit/math/prim", true)
                            runTests("test/unit/math/memory", false)
                        }
                    }
                    post { always { retry(3) { deleteDir() } } }
                }
                stage('OpenCL GPU tests') {
                    agent {
                        docker {
                            image 'stanorg/ci:gpu-cpp17'
                            label 'v100'
                            args '--gpus 1'
                        }
                    }
                    steps {
                        script {
                            unstash 'MathSetup'
                            sh """
                                echo CXX=${CLANG_CXX} -Werror > make/local
                                echo STAN_OPENCL=true >> make/local
                                echo OPENCL_PLATFORM_ID=${OPENCL_PLATFORM_ID_GPU} >> make/local
                                echo OPENCL_DEVICE_ID=${OPENCL_DEVICE_ID_GPU} >> make/local
                            """
                            if (!(params.optimizeUnitTests || isBranch('develop') || isBranch('master'))) {
                                sh "echo O=1 >> make/local"
                            }
                            runTests("test/unit/math/opencl", false) // TODO(bward): try to enable
                            runTests("test/unit/multiple_translation_units_test.cpp")
                        }
                    }
                }
            }
        }
        stage('Always-run tests') {
            when {
                expression {
                    !skipRemainingStages
                }
            }
            failFast true
            parallel {
                stage('MPI tests') {
                    agent {
                        docker {
                            image 'stanorg/ci:gpu-cpp17'
                            label 'linux'
                        }
                    }
                    steps {
                        unstash 'MathSetup'
                        sh """
                            echo CXX=${MPICXX} > make/local
                            echo CXX_TYPE=gcc >> make/local
                            echo STAN_MPI=true >> make/local
                        """
                        runTests("test/unit/math/prim/functor")
                        runTests("test/unit/math/rev/functor")
                    }
                    post { always { retry(3) { deleteDir() } } }
                }
                stage('Expressions test') {
                    agent {
                        docker {
                            image 'stanorg/ci:gpu-cpp17'
                            label 'linux'
                        }
                    }
                    steps {
                        unstash 'MathSetup'
                        script {
                            sh "echo O=0 > make/local"
                            sh "echo CXX=${CLANG_CXX} -Werror >> make/local"
                            sh "python ./test/code_generator_test.py"
                            sh "python ./test/signature_parser_test.py"
                            sh "python ./test/statement_types_test.py"
                            sh "python ./test/varmat_compatibility_summary_test.py"
                            sh "python ./test/varmat_compatibility_test.py"
                            withEnv(['PATH+TBB=./lib/tbb']) {
                                sh "python ./test/expressions/test_expression_testing_framework.py"
                                try { sh "./runTests.py -j${PARALLEL} test/expressions" }
                                finally { junit 'test/**/*.xml' }
                            }
                            sh "make clean-all"
                            sh "echo STAN_THREADS=true >> make/local"
                            withEnv(['PATH+TBB=./lib/tbb']) {
                                try {
                                    sh "./runTests.py -j${PARALLEL} test/expressions --only-functions reduce_sum map_rect"
				                }
                                finally { junit 'test/**/*.xml' }
                            }
                        }
                    }
                    post { always { deleteDir() } }
                }

                stage('Threading tests') {
                    agent {
                        docker {
                            image 'stanorg/ci:gpu-cpp17'
                            label 'linux'
                        }
                    }
                    steps {
                        script {
                            unstash 'MathSetup'
                            sh "echo CXX=${CLANG_CXX} -Werror > make/local"
                            sh "echo STAN_THREADS=true >> make/local"
                            sh "export STAN_NUM_THREADS=4"
                            if (isBranch('develop') || isBranch('master')) {
                                runTests("test/unit")
                                sh "find . -name *_test.xml | xargs rm"
                            } else {
                                runTests("test/unit -f thread")
                                sh "find . -name *_test.xml | xargs rm"
                                runTests("test/unit -f map_rect")
                                sh "find . -name *_test.xml | xargs rm"
                                runTests("test/unit -f reduce_sum")
                            }
                        }
                    }
                    post { always { retry(3) { deleteDir() } } }
                }
            }
        }

        stage ('Discover changed distribution tests') {
            when {
                expression {
                    !skipRemainingStages
                }
            }
            agent {
                docker {
                    image 'stanorg/ci:gpu-cpp17'
                    label 'linux'
                }
            }
            steps {
                script {
                    retry(3) { checkout scm }
                    if (params.runAllDistributions || isBranch('develop') || isBranch('master')) {
                        changedDistributionTests = sh(script:"python3 test/prob/getDependencies.py --pretend-all", returnStdout:true).trim().readLines()
                    } else {
                        changedDistributionTests = sh(script:"python3 test/prob/getDependencies.py", returnStdout:true).trim().readLines()
                    }
                }
            }
        }

        stage ('Distribution tests') {
            when {
                allOf {
                    expression {
                        !skipRemainingStages
                    }
                    expression {
                        !changedDistributionTests.isEmpty()
                    }
                }
            }
            agent { label 'linux && docker' }
            steps {
                script {
                    def tests = [:]
                    for (f in changedDistributionTests.collate(24)) {
                        def names = f.join(" ")
                        tests["Distribution Tests: ${names}"] = { node ("linux && docker") {
                            deleteDir()
                            docker.image('stanorg/ci:gpu-cpp17').inside {
                                catchError {
                                    unstash 'MathSetup'
                                    sh """
                                        echo CXX=${CLANG_CXX} > make/local
                                        echo O=0 >> make/local
                                        echo N_TESTS=${N_TESTS} >> make/local
                                        """
                                    script {
                                        if (params.withRowVector || isBranch('develop') || isBranch('master')) {
                                            sh "echo CXXFLAGS+=-DSTAN_TEST_ROW_VECTORS >> make/local"
                                            sh "echo CXXFLAGS+=-DSTAN_PROB_TEST_ALL >> make/local"
                                        }
                                    }
                                    sh "./runTests.py -j${PARALLEL} ${names}"
                                }
                                deleteDir()
                            }
                        } }
                    }
                    tests.failFast = true
                    parallel tests
                }
            }
            post {
                always {
                    retry(3) { deleteDir() }
                }
            }
        }

        stage('Upstream tests') {
            agent { label 'linux' }
            when {
                allOf {
                    expression {
                        env.BRANCH_NAME ==~ /PR-\d+/
                    }
                    expression {
                        !skipRemainingStages
                    }
                }
            }
            steps {
                build(job: "Stan/Stan/${stan_pr()}",
                        parameters: [string(name: 'math_pr', value: env.BRANCH_NAME),
                                    string(name: 'cmdstan_pr', value: cmdstan_pr())])
            }
        }

        stage('Upload doxygen') {
            agent {
                docker {
                    image 'stanorg/ci:gpu-cpp17'
                    label 'linux'
                }
            }
            when { branch 'develop'}
            steps {
                retry(3) { checkout scm }
                withCredentials([usernamePassword(credentialsId: 'a630aebc-6861-4e69-b497-fd7f496ec46b',
                                                  usernameVariable: 'GIT_USERNAME', passwordVariable: 'GIT_PASSWORD')]) {
                    sh """#!/bin/bash
                        set -x
                        make doxygen
                        git config user.email "mc.stanislaw@gmail.com"
                        git config user.name "Stan Jenkins"
                        git checkout --detach
                        git branch -D gh-pages
                        git push https://${GIT_USERNAME}:${GIT_PASSWORD}@github.com/stan-dev/math.git :gh-pages
                        git checkout --orphan gh-pages
                        git add -f doc
                        git commit --author='Stan BuildBot <mc.stanislaw@gmail.com>' -m "auto generated docs from Jenkins"
                        git subtree push --prefix doc/api/html https://${GIT_USERNAME}:${GIT_PASSWORD}@github.com/stan-dev/math.git gh-pages
                        """
                }
            }
            post { always { deleteDir() } }
        }

    }

    post {
        always {
            node("linux") {
                recordIssues enabledForFailure: false, tool: clang()
            }
        }
        success {
            script {
                utils.updateUpstream(env, 'stan')
                utils.mailBuildResults("SUCCESSFUL")
            }
        }
        unstable { script { utils.mailBuildResults("UNSTABLE", alsoNotify()) } }
        failure { script { utils.mailBuildResults("FAILURE", alsoNotify()) } }
    }
}

#!/usr/bin/env groovy
@Library('StanUtils')
import org.stan.Utils

def runTests(String testPath) {
    sh "./runTests.py -j${env.PARALLEL} ${testPath} --make-only"
    try { sh "./runTests.py -j${env.PARALLEL} ${testPath}" }
    finally { junit 'test/**/*.xml' }
}

def runTestsWin(String testPath) {
    withEnv(['PATH+TBB=./lib/tbb']) {
       bat "echo $PATH"
       bat "runTests.py -j40 ${testPath} --make-only"
       try { bat "runTests.py -j${env.PARALLEL} ${testPath}" }
       finally { junit 'test/**/*.xml' }
    }
}

def deleteDirWin() {
    bat "attrib -r -s /s /d"
    deleteDir()
}

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
        string(defaultValue: '', name: 'cmdstan_pr',
          description: 'PR to test CmdStan upstream against e.g. PR-630')
        string(defaultValue: '', name: 'stan_pr',
          description: 'PR to test Stan upstream against e.g. PR-630')
        booleanParam(defaultValue: false, description:
        'Run additional distribution tests on RowVectors (takes 5x as long)',
        name: 'withRowVector')
    }
    options {
        skipDefaultCheckout()
        preserveStashes(buildCount: 7)
    }
    environment {
        STAN_NUM_THREADS = '4'
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
        stage('Always-run tests part 2') {
            parallel {
                stage('Windows Headers & Unit') {
                    agent { label 'windows' }
                    steps {
                        deleteDirWin()
                        unstash 'MathSetup'
                        bat "mingw32-make -j${env.PARALLEL} test-headers"
                        runTestsWin("test/unit")
                    }
                }
                stage('Windows Threading') {
                    agent { label 'windows' }
                    steps {
                        deleteDirWin()
                        unstash 'MathSetup'
                        bat "echo CXX=${env.CXX} -Werror > make/local"
                        bat "echo CXXFLAGS+=-DSTAN_THREADS >> make/local"
                        runTestsWin("test/unit -f thread")
                        runTestsWin("test/unit -f map_rect")
                    }
                }
            }
        }
    }
    post {
        always {
            node("osx || linux") {
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

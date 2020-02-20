#!/usr/bin/env groovy
@Library('StanUtils')
import org.stan.Utils

def runTests(String testPath) {
    sh "./runTests.py -j${env.PARALLEL} ${testPath} --make-only"
    try { sh "./runTests.py -j${env.PARALLEL} ${testPath}" }
    finally { junit 'test/**/*.xml' }
}

def runTestsWin(String testPath, boolean buildLibs = true) {
    withEnv(['PATH+TBB=./lib/tbb']) {
       bat "echo $PATH"
       if (buildLibs){
           bat "mingw32-make.exe -f make/standalone math-libs"
       }
       bat "runTests.py -j${env.PARALLEL} ${testPath} --make-only"
       try { bat "runTests.py -j${env.PARALLEL} ${testPath}" }
       finally { junit 'test/**/*.xml' }
    }
}

def deleteDirWin() {
    bat "attrib -r -s /s /d"
    deleteDir()
}

def sourceCodePaths(){
    // These paths will be passed to git diff
    // If there are changes to them, CI/CD will continue else skip
    def paths = ['stan', 'make', 'lib', 'test', 'runTests.py', 'runChecks.py', 'makefile', 'Jenkinsfile', '.clang-format']
    def bashArray = ""

    for(path in paths){
        bashArray += path + (path != paths[paths.size() - 1] ? " " : "")
    }

    return bashArray
}

def skipRemainingStages = true

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
        scPaths = sourceCodePaths()
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
            agent any
            steps {
                sh "printenv"
                deleteDir()
                retry(3) { checkout scm }
                withCredentials([usernamePassword(credentialsId: 'a630aebc-6861-4e69-b497-fd7f496ec46b',
                    usernameVariable: 'GIT_USERNAME', passwordVariable: 'GIT_PASSWORD')]) {
                    sh """#!/bin/bash
                        set -x
                        git checkout -b ${branchName()}
                        clang-format --version
                        find stan test -name '*.hpp' -o -name '*.cpp' | xargs -n20 -P${env.PARALLEL} clang-format -i
                        if [[ `git diff` != "" ]]; then
                            git config --global user.email "mc.stanislaw@gmail.com"
                            git config --global user.name "Stan Jenkins"
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
            agent any
            steps {
                script {
                    deleteDir()
                    retry(3) { checkout scm }
                    sh "git clean -xffd"
                    stash 'MathSetup'
                    sh "echo CXX=${env.CXX} -Werror > make/local"
                    sh "echo BOOST_PARALLEL_JOBS=${env.PARALLEL} >> make/local"
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
            agent { label 'linux' }
            steps {
                script {         

                    retry(3) { checkout scm }
                    sh 'git clean -xffd'

                    def commitHash = sh(script: "git rev-parse HEAD | tr '\\n' ' '", returnStdout: true)
                    def changeTarget = ""

                    if (env.CHANGE_TARGET) {
                        println "This build is a PR, checking out target branch to compare changes."
                        changeTarget = env.CHANGE_TARGET
                        sh(script: "git pull && git checkout ${changeTarget}", returnStdout: false)
                    }
                    else{
                        println "This build is not PR, checking out current branch and extract HEAD^1 commit to compare changes or develop when downstream_tests."
                        if (env.BRANCH_NAME == "downstream_tests"){
                            sh(script: "git checkout develop && git pull", returnStdout: false)
                            changeTarget = sh(script: "git rev-parse HEAD^1 | tr '\\n' ' '", returnStdout: true)
                            sh(script: "git checkout ${commitHash}", returnStdout: false)
                        }
                        else{
                            sh(script: "git pull && git checkout ${env.BRANCH_NAME}", returnStdout: false)
                            changeTarget = sh(script: "git rev-parse HEAD^1 | tr '\\n' ' '", returnStdout: true)
                        }
                    }

                    println "Comparing differences between current ${commitHash} and target ${changeTarget}"

                    def bashScript = """
                        for i in ${env.scPaths};
                        do
                            git diff ${commitHash} ${changeTarget} -- \$i
                        done
                    """

                    def differences = sh(script: bashScript, returnStdout: true)

                    println differences

                    if (differences?.trim()) {
                        println "There are differences in the source code, CI/CD will run."
                        skipRemainingStages = false
                    }
                    else{
                        println "There aren't any differences in the source code, CI/CD will not run."
                        skipRemainingStages = true
                    }
                }
            }
        }
        stage('Headers checks') {
            when {
                expression {
                    !skipRemainingStages
                }
            }
            parallel {
              stage('Headers check') {
                agent any
                steps {
                    deleteDir()
                    unstash 'MathSetup'
                    sh "echo CXX=${env.CXX} -Werror > make/local"
                    sh "make -j${env.PARALLEL} test-headers"
                }
                post { always { deleteDir() } }
              }
              stage('Headers check with OpenCL') {
                agent { label "gpu" }
                steps {
                    deleteDir()
                    unstash 'MathSetup'
                    sh "echo CXX=${env.CXX} -Werror > make/local"
                    sh "echo STAN_OPENCL=true>> make/local"
                    sh "echo OPENCL_PLATFORM_ID=0>> make/local"
                    sh "echo OPENCL_DEVICE_ID=${OPENCL_DEVICE_ID}>> make/local"
                    sh "make -j${env.PARALLEL} test-headers"
                }
                post { always { deleteDir() } }
              }
           }
        }
        stage('Always-run tests part 1') {
            when {
                expression {
                    !skipRemainingStages
                }
            }
            parallel {
                stage('Linux Unit with MPI') {
                    agent { label 'linux && mpi' }
                    steps {
                        deleteDir()
                        unstash 'MathSetup'
                        sh "echo CXX=${MPICXX} >> make/local"
                        sh "echo CXX_TYPE=gcc >> make/local"                        
                        sh "echo STAN_MPI=true >> make/local"
                        runTests("test/unit")
                    }
                    post { always { retry(3) { deleteDir() } } }
                }
                stage('Full unit with GPU') {
                    agent { label "gpu" }
                    steps {
                        deleteDir()
                        unstash 'MathSetup'
                        sh "echo CXX=${env.CXX} -Werror > make/local"
                        sh "echo STAN_OPENCL=true>> make/local"
                        sh "echo OPENCL_PLATFORM_ID=0>> make/local"
                        sh "echo OPENCL_DEVICE_ID=${OPENCL_DEVICE_ID}>> make/local"
                        runTests("test/unit")
                    }
                    post { always { retry(3) { deleteDir() } } }
                }
            }
        }
        stage('Always-run tests part 2') {
            when {
                expression {
                    !skipRemainingStages
                }
            }
            parallel {
                stage('Distribution tests') {
                    agent { label "distribution-tests" }
                    steps {
                        deleteDir()
                        unstash 'MathSetup'
                        sh """
                            echo CXX=${env.CXX} > make/local
                            echo O=0 >> make/local
                            echo N_TESTS=${env.N_TESTS} >> make/local
                            """
                        script {
                            if (params.withRowVector || isBranch('develop') || isBranch('master')) {
                                sh "echo CXXFLAGS+=-DSTAN_TEST_ROW_VECTORS >> make/local"
                            }
                        }
                        sh "./runTests.py -j${env.PARALLEL} test/prob > dist.log 2>&1"
                    }
                    post {
                        always {
                            script { zip zipFile: "dist.log.zip", archive: true, glob: 'dist.log' }
                            retry(3) { deleteDir() }
                        }
                        failure {
                            echo "Distribution tests failed. Check out dist.log.zip artifact for test logs."
                            }
                    }
                }
                stage('Threading tests') {
                    agent any
                    steps {
                        deleteDir()
                        unstash 'MathSetup'
                        sh "echo CXX=${env.CXX} -Werror > make/local"
                        sh "echo CPPFLAGS+=-DSTAN_THREADS >> make/local"
                        runTests("test/unit -f thread")
                        sh "find . -name *_test.xml | xargs rm"
                        runTests("test/unit -f map_rect")
                    }
                    post { always { retry(3) { deleteDir() } } }
                }
                stage('Windows Headers & Unit') {
                    agent { label 'windows' }
                    steps {
                        deleteDirWin()
                        unstash 'MathSetup'
                        bat "mingw32-make.exe -f make/standalone math-libs"
                        bat "mingw32-make -j${env.PARALLEL} test-headers"
                        runTestsWin("test/unit", false)
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
        stage('Additional merge tests') {
            when { 
                allOf {
                    anyOf {
                        branch 'develop'
                        branch 'master'
                    }
                    expression {
                        !skipRemainingStages
                    }
                }
            }
            parallel {
                stage('Linux Unit with Threading') {
                    agent { label 'linux' }
                    steps {
                        deleteDir()
                        unstash 'MathSetup'
                        sh "echo CXX=${GCC} >> make/local"
                        sh "echo CXXFLAGS=-DSTAN_THREADS >> make/local"
                        runTests("test/unit")
                    }
                    post { always { retry(3) { deleteDir() } } }
                }
                stage('Mac Unit with Threading') {
                    agent  { label 'osx' }
                    steps {
                        deleteDir()
                        unstash 'MathSetup'
                        sh "echo CC=${env.CXX} -Werror > make/local"
                        sh "echo CXXFLAGS+=-DSTAN_THREADS >> make/local"
                        runTests("test/unit")
                    }
                    post { always { retry(3) { deleteDir() } } }
                }
            }
        }
        stage('Upstream tests') {
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
                build(job: "Stan/${stan_pr()}",
                        parameters: [string(name: 'math_pr', value: env.BRANCH_NAME),
                                    string(name: 'cmdstan_pr', value: cmdstan_pr())])
            }
        }
        stage('Upload doxygen') {
            agent any
            when { branch 'develop'}
            steps {
                deleteDir()
                retry(3) { checkout scm }
                withCredentials([usernamePassword(credentialsId: 'a630aebc-6861-4e69-b497-fd7f496ec46b',
                                                  usernameVariable: 'GIT_USERNAME', passwordVariable: 'GIT_PASSWORD')]) {
                    sh """#!/bin/bash
                        set -x
                        make doxygen
                        git config --global user.email "mc.stanislaw@gmail.com"
                        git config --global user.name "Stan Jenkins"
                        git checkout --detach
                        git branch -D gh-pages
                        git push https://${GIT_USERNAME}:${GIT_PASSWORD}@github.com/stan-dev/math.git :gh-pages
                        git checkout --orphan gh-pages
                        git add -f doc
                        git commit -m "auto generated docs from Jenkins"
                        git subtree push --prefix doc/api/html https://${GIT_USERNAME}:${GIT_PASSWORD}@github.com/stan-dev/math.git gh-pages
                        """
                }
            }
            post { always { deleteDir() } }
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

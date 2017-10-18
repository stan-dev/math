def setupCC(Boolean failOnError = true) {
    errorStr = failOnError ? "-Werror " : ""
    "echo CC=${env.CXX} ${errorStr}> make/local"
}

def setup(Boolean failOnError = true) {
    sh """
        git clean -xffd
        ${setupCC(failOnError)}
    """
}

def mailBuildResults(String label, additionalEmails='') {
    emailext (
        subject: "[StanJenkins] ${label}: Job '${env.JOB_NAME} [${env.BUILD_NUMBER}]'",
        body: """${label}: Job '${env.JOB_NAME} [${env.BUILD_NUMBER}]': Check console output at ${env.BUILD_URL}""",
        recipientProviders: [[$class: 'RequesterRecipientProvider']],
        to: "${env.CHANGE_AUTHOR_EMAIL}, ${additionalEmails}"
    )
}

def runTests(String testPath) {
    sh "./runTests.py -j${env.PARALLEL} ${testPath} || echo ${testPath} failed"
}

def updateUpstream(String upstreamRepo) {
    if (env.BRANCH_NAME == 'develop') {
        node('master') {
            retry(3) {
                checkout scm
            }
            sh """
                curl -O https://raw.githubusercontent.com/stan-dev/ci-scripts/master/jenkins/create-${upstreamRepo}-pull-request.sh
                sh create-${upstreamRepo}-pull-request.sh
            """
            deleteDir()
        }
    }
}

pipeline {
    agent none
    parameters {
        string(defaultValue: 'downstream tests', name: 'cmdstan_pr',
          description: 'PR to test CmdStan upstream against e.g. PR-630')
        string(defaultValue: 'downstream tests', name: 'stan_pr',
          description: 'PR to test Stan upstream against e.g. PR-630')
    }
    options { skipDefaultCheckout() }
    stages {
        stage('Linting & Doc checks') {
            agent any
            steps {
                script {
                    setup(false)
                    stash 'MathSetup'
                    sh setupCC()
                    parallel(
                        CppLint: { sh "make cpplint" },
                        dependencies: { sh 'make test-math-dependencies' } ,
                        documentation: { sh 'make doxygen' },
                        failFast: true
                    )
                }
                post { always { deleteDir() } }
            }
        }
        stage('Tests') {
            failFast true
            parallel {
                stage('Headers') {
                    agent any
                    steps { 
                        unstash 'MathSetup'
                        sh setupCC()
                        sh "make -j${env.PARALLEL} test-headers"
                    }
                    post { always { deleteDir() } }
                }
                stage('Unit') {
                    agent any
                    steps {
                        unstash 'MathSetup'
                        sh setupCC()
                        runTests("test/unit")
                        retry(2) { junit 'test/**/*.xml' }
                    }
                    post { always { deleteDir() } }
                }
                stage('CmdStan Upstream Tests') {
                    when { expression { env.BRANCH_NAME ==~ /PR-\d+/ } }
                    steps {
                        build(job: "CmdStan/${params.cmdstan_pr}",
                                    parameters: [string(name: 'math_pr', value: env.BRANCH_NAME)])
                    }
                }
                stage('Stan Upstream Tests') {
                    when { expression { env.BRANCH_NAME ==~ /PR-\d+/ } }
                    steps {
                        build(job: "Stan Pipeline/${params.stan_pr}",
                                    parameters: [string(name: 'math_pr', value: env.BRANCH_NAME)])
                    }
                }
                stage('Distribution tests') {
                    agent { label "distribution-tests" }
                    // XXX Add conditional back in so we don't run this if we haven't
                    // changed code or makefiles
                    steps { 
                        unstash 'MathSetup'
                        sh setupCC(false)
                        runTests("test/prob")
                        retry(2) { junit 'test/**/*.xml' }
                    }
                    post { always { deleteDir() } }
                }
            }
        }
    }
    post {
        always {
            node('master') {
                warnings consoleParsers: [[parserName: 'GNU C Compiler 4 (gcc)']], canRunOnFailed: true
                warnings consoleParsers: [[parserName: 'Clang (LLVM based)']], canRunOnFailed: true
                warnings consoleParsers: [[parserName: 'CppLint']], canRunOnFailed: true
                warnings consoleParsers: [[parserName: 'math-dependencies']], canRunOnFailed: true
            }
        }
        success {
            updateUpstream('stan')
            mailBuildResults("SUCCESSFUL")
        }
        unstable { mailBuildResults("UNSTABLE", "stan-buildbot@googlegroups.com") }
        failure { mailBuildResults("FAILURE", "stan-buildbot@googlegroups.com") }
    }
}

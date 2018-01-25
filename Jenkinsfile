@Library('StanUtils')
import org.stan.Utils

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
    sh "./runTests.py -j${env.PARALLEL} ${testPath}"
}

def utils = new org.stan.Utils()

def isBranch(String b) { env.BRANCH_NAME == b }

def alsoNotify() {
    if (isBranch('master') || isBranch('develop')) {
        "stan-buildbot@googlegroups.com"
    } else ""
}

pipeline {
    agent none
    parameters {
        string(defaultValue: 'downstream tests', name: 'cmdstan_pr',
          description: 'PR to test CmdStan upstream against e.g. PR-630')
        string(defaultValue: 'downstream tests', name: 'stan_pr',
          description: 'PR to test Stan upstream against e.g. PR-630')
        booleanParam(defaultValue: false, description:
        'Run additional distribution tests on RowVectors (takes 5x as long)',
        name: 'withRowVector')
    }
    options { skipDefaultCheckout() }
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
        stage('Lint, doc, header, unit tests') {
            agent any
            steps { script {
                stage('Git checkout') {
                    sh "printenv"
                    retry(3) { checkout scm }
                    setup(false)
                    stash 'MathSetup'
                    sh setupCC()
                }
                stage('CppLint') {
                    sh "find test/unit stan -name '*.hpp' -or -name '*.cpp' " +
                        "| xargs -n20 -P${env.PARALLEL}" +
                        " python lib/cpplint_4.45/cpplint.py --output=vs7 --counting=detailed --root=stan --extension=hpp,cpp --filter=-runtime/indentation_namespace,-build/c++11,-readability/namespace,-legal/copyright,-whitespace/indent,-runtime/reference"
                }
                stage('Dependencies') {
                    sh "make test-math-dependencies"
                }
                stage('Documentation') {
                    sh 'make doxygen'
                }
                stage('Headers tests') {
                    sh "make -j${env.PARALLEL} test-headers"
                }
                stage('Unit tests') {
                    runTests("test/unit")
                }
            } }
            post {
                always {
                    script {
                        if (env.CXX.contains("clang")) {
                            warnings consoleParsers: [[parserName: 'Clang (LLVM based)']], canRunOnFailed: true
                        } else {
                            warnings consoleParsers: [[parserName: 'GNU C Compiler 4 (gcc)']], canRunOnFailed: true
                        }
                    }
                    warnings consoleParsers: [[parserName: 'math-dependencies']], canRunOnFailed: true
                    warnings consoleParsers: [[parserName: 'CppLint']], canRunOnFailed: true
                    junit 'test/**/*.xml'
                    deleteDir()
                }
            }
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
                build(job: "Stan/${params.stan_pr}",
                            parameters: [string(name: 'math_pr', value: env.BRANCH_NAME)])
            }
        }
        stage('Distribution tests') {
            agent { label "distribution-tests" }
            steps {
                unstash 'MathSetup'
                sh """
                    ${setupCC(false)}
                    echo 'O=0' >> make/local
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
    }
    post {
        success {
            script { utils.updateUpstream('stan') }
            mailBuildResults("SUCCESSFUL")
        }
        unstable { mailBuildResults("UNSTABLE", alsoNotify()) }
        failure {
            mailBuildResults("FAILURE", alsoNotify())
        }
    }
}

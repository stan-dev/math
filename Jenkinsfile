def clean() {
    sh """
        git clean -xffd
        echo 'CC=${env.CXX}' > make/local
    """
}

def mailDevs(String label) {
    emailext (
        subject: "[StanJenkins] ${label}: Job '${env.JOB_NAME} [${env.BUILD_NUMBER}]'",
        body: """<p>${label}: Job '${env.JOB_NAME} [${env.BUILD_NUMBER}]':</p>
            <p>Check console output at &QUOT;<a href='${env.BUILD_URL}'>${env.JOB_NAME} [${env.BUILD_NUMBER}]</a>&QUOT;</p>""",
        recipientProviders: [[$class: 'DevelopersRecipientProvider']]
    )
}

def runTests(String testPath) {
    sh "./runTests.py -j${env.PARALLEL} ${testPath} || echo ${testPath} failed"
}

pipeline {
    agent none
    parameters {
        string(defaultValue: 'downstream tests', name: 'cmdstan_pr',
          description: 'PR to test CmdStan upstream against e.g. PR-630')
        string(defaultValue: 'downstream tests', name: 'stan_pr',
          description: 'PR to test Stan upstream against e.g. PR-630')
    }
    stages {
        stage('Linting & Doc checks') {
            agent any
            steps {
                script {
                    clean()
                    sh "echo 'CXXFLAGS += -Werror' >> make/local"
                    parallel(
                        CppLint: { sh "make cpplint" },
                        dependencies: { sh 'make test-math-dependencies' } ,
                        documentation: { sh 'make doxygen' },
                        failFast: true
                    )
                }
            }
        }
        stage('Tests') {
            failFast true
            parallel {
                stage('Headers') {
                    agent any
                    steps { 
                        clean()
                        sh "echo 'CXXFLAGS += -Werror' >> make/local"
                        sh "make -j${env.PARALLEL} test-headers"
                    }
                }
                stage('Unit') {
                    agent any
                    steps {
                        clean()
                        sh "echo 'CXXFLAGS += -Werror' >> make/local"
                        runTests("test/unit")
                        junit 'test/**/*.xml'
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
                        build(job: "Stan Pipeline/${params.stan_pr}",
                                    parameters: [string(name: 'math_pr', value: env.BRANCH_NAME)])
                    }
                }
                stage('Distribution tests') {
                    agent { label "distribution-tests" }
                    // XXX Add conditional back in so we don't run this if we haven't
                    // changed code or makefiles
                    steps { 
                        clean()
                        runTests("test/prob")
                        junit 'test/**/*.xml'
                    }
                    post { always { cleanWs() } }
                }
            }
        }
        stage('Update Stan Upstream') {
           agent none
            when { branch "develop" }
            steps {
                sh "curl -O https://raw.githubusercontent.com/stan-dev/ci-scripts/master/jenkins/create-stan-pull-request.sh"
                sh "sh create-stan-pull-request.sh"
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
        success { mailDevs("SUCCESSFUL") }
        failure { mailDevs("FAILURE") }
    }
}

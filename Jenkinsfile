def clean() {
    sh """
        git clean -xffd
        echo 'CC=${env.CXX}' > make/local
    """
}

pipeline {
    agent { label 'master' }
    options {
        disableConcurrentBuilds()
    }
    parameters {
        string(defaultValue: 'develop', name: 'cmdstan_pr',
          description: 'PR to test CmdStan against e.g. PR-630')
    }
    stages {
        stage('Clean & Setup') {
            steps {
                clean()
                sh "echo 'CXXFLAGS += -Werror' >> make/local"
            }
        }
        stage('Linting & Doc checks') {
            steps {
                parallel(
                    CppLint: { sh "make cpplint" },
                    dependencies: { sh 'make test-math-dependencies' } ,
                    documentation: { sh 'make doxygen' },
                    failFast: true
                )
            }
        }
        stage('Headers') {
            steps { sh "make -j${env.PARALLEL} test-headers" }
        }
        stage('Unit') {
            agent any
            steps {
                clean()
                sh "echo 'CXXFLAGS += -Werror' >> make/local"
                sh "./runTests.py -j${env.PARALLEL} test/unit"
            }
        }
        stage('Upstream Tests') {
            agent none
            when { 
                allOf {
                    not { branch 'master' }
                    not { branch 'develop' }
                }
            }
            steps {
                parallel(
                    //stan: { build(job: 'Stan Pipeline',
                    //            parameters: [string(name: 'math_pr', value: env.BRANCH_NAME)])},
                    cmdStan: { build(job: "CmdStan Pipeline/${params.cmdstan_pr}",
                                parameters: [string(name: 'math_pr', value: env.BRANCH_NAME)])},
                    failFast: true
                )
            }
        }
        stage('Distribution tests') {
            agent { label "distribution-tests" }
            // XXX Add conditional back in so we don't run this if we haven't
            // changed code or makefiles
            steps { 
                clean()
                sh "./runTests.py -j${env.PARALLEL} test/prob"
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
            warnings consoleParsers: [[parserName: 'GNU C Compiler 4 (gcc)']], canRunOnFailed: true
            warnings consoleParsers: [[parserName: 'Clang (LLVM based)']], canRunOnFailed: true
            warnings consoleParsers: [[parserName: 'CppLint']], canRunOnFailed: true
            warnings consoleParsers: [[parserName: 'math-dependencies']], canRunOnFailed: true
        }
    }
}

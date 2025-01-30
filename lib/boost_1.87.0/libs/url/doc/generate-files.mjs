#!/usr/bin/env node

//
// Copyright (c) 2023 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

'use strict';
import {execSync} from 'child_process'
import fs from 'fs'
import path from 'path'
import os from 'os'
import request from 'sync-request'

/*
    Helper functions
 */
function mkdirp(dir) {
    if (!fs.existsSync(dir)) {
        mkdirp(path.dirname(dir))
        fs.mkdirSync(dir)
    }
}

function findExecutable(executableName) {
    if (Array.isArray(executableName)) {
        for (const name of executableName) {
            const result = findExecutable(name)
            if (result) {
                return result
            }
        }
        return undefined
    }

    const isWin = process.platform === 'win32';
    const pathDirs = process.env.PATH.split(isWin ? ';' : ':');
    const extensions = isWin ? ['.exe', '.bat', '.cmd'] : [''];

    function isExecutable(filePath) {
        try {
            if (!isWin) {
                fs.accessSync(filePath, fs.constants.X_OK);
            }
            return true;
        } catch (error) {
            return false;
        }
    }

    // Try to find the exact executable first
    for (const dir of pathDirs) {
        for (const ext of extensions) {
            const fullPath = path.join(dir, executableName + ext);
            if (fs.existsSync(fullPath) && isExecutable(fullPath)) {
                return fullPath;
            }
        }
    }

    function escapeRegExp(string) {
        // Escape special characters for use in regex
        return string.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
    }

    // If the exact executable is not found, search for versioned executables
    const versionedExecutables = [];
    const escapedExecutableName = escapeRegExp(executableName);
    const versionRegex = new RegExp(`${escapedExecutableName}-(\\d+)$`);

    for (const dir of pathDirs) {
        try {
            const files = fs.readdirSync(dir);
            for (const file of files) {
                if (!extensions.some(ext => file.endsWith(ext))) {
                    continue
                }
                const fullPath = path.join(dir, file);
                if (!isExecutable(fullPath)) {
                    continue
                }
                const ext = path.extname(file);
                const basename = path.basename(file, ext);
                const match = basename.match(versionRegex);
                if (match) {
                    versionedExecutables.push({
                        path: fullPath,
                        version: parseInt(match[1], 10)
                    });
                }
            }
        } catch (error) {
            // Ignore errors from reading directories
        }
    }

    if (versionedExecutables.length > 0) {
        versionedExecutables.sort((a, b) => b.version - a.version);
        return versionedExecutables[0].path;
    }

    return undefined;
}

function mkTmpDir() {
    let tempDir = undefined
    try {
        tempDir = fs.mkdtempSync(os.tmpdir())
    } catch (error) {
        console.error(`Could not create temporary directory with fs.mkdtempSync(os.tmpdir()): ${error}.`)
        try {
            console.error(`Trying to create temporary directory in fs.mkdtempSync('~/tmp')...`)
            tempDir = fs.mkdtempSync(path.join(os.homedir(), 'tmp'));
        } catch (error) {
            console.error(`Could not create temporary directory with fs.mkdtempSync('~/tmp'): ${error}.`)
            try {
                console.error(`Trying to create temporary directory in fs.mkdtempSync(path.join(cwd, 'tmp'))...`)
                tempDir = fs.mkdtempSync(path.join(cwd, 'tmp'));
                console.log(`Created temporary directory ${tempDir}`)
            } catch (error) {
                console.error(`Could not create temporary directory: ${error}`)
            }
        }
    }
    return tempDir
}

function downloadAndDecompress(downloadUrl, downloadPath, extractPath) {
    try {
        console.log(`Downloading ${downloadUrl} to ${downloadPath}...`);
        const response = request('GET', downloadUrl)
        fs.writeFileSync(downloadPath, response.getBody());
        console.log(`File ${downloadUrl} downloaded successfully to ${downloadPath}.`);
        if (path.extname(downloadPath) === '.7z') {
            const cmd = `7z x ${downloadPath} -o${extractPath}`
            console.log(`Decompressing with command: ${cmd}`)
            execSync(cmd, {stdio: 'inherit'})
        } else if (/\.tar\.gz$/.test(downloadPath)) {
            mkdirp(extractPath)
            const cmd = `tar -vxzf ${downloadPath} -C ${extractPath}`
            console.log(`Decompressing with command: ${cmd}`)
            execSync(cmd, {stdio: 'inherit'})
        }
        console.log(`File decompressed successfully to ${extractPath}.`);
    } catch (error) {
        console.error('Error:', error.message);
    }
}

/*
    Find C++ compilers
    Clang++ is preferred over g++ and cl for MrDocs
 */
const cwd = process.cwd()
let cxxCompiler = findExecutable(['clang++', 'g++', 'cl']) || process.env.CXX_COMPILER || process.env.CXX
if (cxxCompiler && process.platform === "win32") {
    // Replace "\" with "/" in CXX_COMPILER
    cxxCompiler = cxxCompiler.replace(/\\/g, '/')
}
if (cxxCompiler === undefined) {
    console.error('Could not find a C++ compiler. Please set the CXX_COMPILER environment variable.')
    process.exit(1)
}
process.env.CMAKE_CXX_COMPILER = cxxCompiler
process.env.CXX = cxxCompiler

const cxxCompilerName = path.basename(cxxCompiler).replace(/\.exe$/, '')
let cCompiler = findExecutable(['clang', 'gcc', 'cl']) || process.env.C_COMPILER || process.env.CC
if (cCompiler && process.platform === "win32") {
    // Replace "\" with "/" in CXX_COMPILER
    cCompiler = cCompiler.replace(/\\/g, '/')
}
if (cCompiler === undefined) {
    console.error('Could not find a C compiler. Please set the C_COMPILER environment variable.')
    process.exit(1)
}
const cCompilerName = path.basename(cCompiler).replace(/\.exe$/, '')
process.env.CMAKE_C_COMPILER = cCompiler
process.env.CC = cCompiler

console.log(`C++ compiler: ${cxxCompilerName} (${cxxCompiler})`)
console.log(`C compiler: ${cCompilerName} (${cCompiler})`)

/*
    Download Boost
 */
function isBoostDir(dir) {
    if (!fs.existsSync(dir)) {
        return false
    }
    console.log(`Checking if ${dir} is a Boost source directory`)
    const boostFiles = ['CMakeLists.txt', 'Jamroot', 'boost-build.jam', 'bootstrap.sh', 'libs']
    for (const file of boostFiles) {
        const filePath = path.join(cwdParentParent, file)
        console.log(`Checking if ${filePath} exists`)
        if (!fs.existsSync(filePath)) {
            return false
        }
    }
    return true
}

// Check if we need to download boost. This is usually the case in CI, but we want to reuse the local boost directory
// in local builds.
const librariesJsonPath = path.join(cwd, 'meta', 'libraries.json')
const librariesJson = JSON.parse(fs.readFileSync(librariesJsonPath, 'utf-8'))
const module = librariesJson['key']
const self = module || path.basename(cwd)
console.log(`self: ${self}`)

// Check if cwdParentParent contains the files "CMakeLists.txt" "Jamroot" "boost-build.jam" "bootstrap.sh" "libs"
let tempDir = undefined
let boostDir = undefined
const cwdParentParent = path.dirname(path.dirname(cwd))
if (isBoostDir(cwdParentParent)) {
    boostDir = cwdParentParent
    console.log(`Found Boost at ${cwdParentParent}`)
} else {
    const branch = execSync('git rev-parse --abbrev-ref HEAD', {encoding: 'utf-8'}).trim()
    const boostBranchRegex = /^boost-1\.\d+\.\d+$/
    const isValidBoostBranch = branch === 'master' || branch === 'develop' || boostBranchRegex.test(branch);
    const boostBranch = isValidBoostBranch ? branch : 'develop'

    // Clone boost in a temporary directory
    if (!tempDir) {
        tempDir = mkTmpDir()
    }
    console.log(`Created temporary directory ${tempDir}`)
    let boostDownloadCmd = `git clone https://github.com/boostorg/boost.git --depth 1 --branch ${boostBranch}`
    execSync(boostDownloadCmd, {cwd: tempDir})
    boostDir = path.join(tempDir, 'boost')
    console.log(`Initializing submodules in ${boostDir}`)
    execSync('git submodule update --init --recursive', {cwd: boostDir})
    console.log(`Cloned Boost to ${boostDir}`)

    // Delete `self` from boost/libs
    const selfDir = path.join(boostDir, 'libs', self)
    console.log(`Deleting ${selfDir}`)
    execSync(`rm -rf ${selfDir}`)

    // Copy contents of cwd to boost/libs/self
    // const selfDirParent = path.join(boostDir, 'libs')
    console.log(`Copying ${cwd} to ${selfDir}`)
    execSync(`cp -r ${cwd} ${selfDir}`)
}

/*
    Find and run MrDocs to generate the documentation

    Install MrDocs with:
        - Linux:
            wget https://github.com/cppalliance/mrdocs/releases/download/develop-release/MrDocs-0.0.1-Linux.tar.gz
            sudo tar -xzf MrDocs-0.0.1-Linux.tar.gz -C /usr/local --strip-components=1
        - Windows (Powershell):
            Invoke-WebRequest -Uri 'https://github.com/cppalliance/mrdocs/releases/download/develop-release/MrDocs-0.0.1-win64.7z' -OutFile 'MrDocs-0.0.1-win64.7z'
            7z x -o"C:\Users\$env:USERNAME\Applications" 'MrDocs-0.0.1-win64.7z'
            ( Adapt the destination path as needed )

 */
// https://github.com/terascope/fetch-github-release
let mrDocsExec = findExecutable('mrdocs')
if (!mrDocsExec) {
    console.log(`Could not find MrDocs. Downloading...`)
    if (!tempDir) {
        tempDir = mkTmpDir()
    }
    const releasesResponse = request('GET', 'https://api.github.com/repos/cppalliance/mrdocs/releases', {
        headers: {
            'User-Agent': 'request'
        }
    })
    const releasesInfo = JSON.parse(releasesResponse.getBody('utf-8'))
    console.log(`Found ${releasesInfo.length} MrDocs releases`)
    let downloadUrl = undefined
    for (const latestRelease of releasesInfo) {
        console.log(`Latest release: ${latestRelease['tag_name']}`)
        const latestAssets = latestRelease['assets'].map(asset => asset['browser_download_url'])
        console.log(`Latest assets: ${latestAssets}`)
        downloadUrl = process.platform === "win32" ? latestAssets.find(asset => asset.endsWith('win64.7z')) : latestAssets.find(asset => asset.endsWith('Linux.tar.gz'))
        if (downloadUrl) {
            break
        }
        console.warn(`Could not find MrDocs binaries in ${latestRelease['tag_name']} release for ${process.platform}`)
    }
    if (!downloadUrl) {
        console.error(`Could not find MrDocs binaries for ${process.platform}`)
        process.exit(1)
    }
    const downloadFilename = path.basename(downloadUrl)
    console.log(`Downloading ${downloadUrl} to ${path.join(tempDir, downloadFilename)}...`)
    downloadAndDecompress(
        downloadUrl,
        path.join(tempDir, downloadFilename),
        path.join(tempDir, 'MrDocs')
    )
    console.log(`Extracted ${downloadFilename} to ${path.join(tempDir, 'MrDocs')}`)
    const downloadExtractDir = downloadFilename.replace(/\.(7z|tar\.gz)$/, '')
    const downloadExtractBinFilename = process.platform === "win32" ? 'mrdocs.exe' : 'mrdocs'
    mrDocsExec = path.join(tempDir, 'MrDocs', downloadExtractDir, 'bin', downloadExtractBinFilename)
    if (!fs.existsSync(mrDocsExec)) {
        console.error(`Could not find MrDocs at ${mrDocsExec}`)
        process.exit(1)
    } else {
        console.log(`Found MrDocs executable at ${mrDocsExec}`)
    }
}
console.log(`Found MrDocs at ${mrDocsExec}`)

// Reference goes to another module so that its relative links work
// Antora does not support relative xrefs: https://gitlab.com/antora/antora/-/issues/428
const buildDirectory = path.join(cwd, 'doc', 'build')
const mrDocsOutputDir = path.join(buildDirectory, 'generated-files', 'modules', 'reference', 'pages')
const mrDocsProjectPath = path.join(boostDir, 'libs', 'url')
const mrDocsConfigPath = path.join(mrDocsProjectPath, 'doc', 'mrdocs.yml')
const mrDocsCmd = `${mrDocsExec} --config="${mrDocsConfigPath}" "${mrDocsProjectPath}" --output="${mrDocsOutputDir}"`
console.log(`Generating documentation with command: ${mrDocsCmd}`)
// Execute mrDocsCmd with execSync, get output and print it to stdout. Also get return code.
// If return code is not 0, print output to stderr and exit with return code.
let mrDocsExitCode = 0
try {
    const mrDocsOutput = execSync(mrDocsCmd)
    console.log(mrDocsOutput.toString())
    console.log(`Generated documentation at ${mrDocsOutputDir}`)
} catch (error) {
    console.error('Failed to run MrDocs')
    console.error(error.stdout.toString())
    mrDocsExitCode = error.status
}
if (tempDir) {
    console.log(`Deleting temporary directory ${tempDir}`)
    execSync(`rm -rf ${tempDir}`)
}
if (mrDocsExitCode !== 0) {
    console.log(`Exiting with code ${mrDocsExitCode}`)
}
process.exit(mrDocsExitCode)

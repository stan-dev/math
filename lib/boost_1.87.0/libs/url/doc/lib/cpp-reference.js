/*
    Copyright (c) 2024 Alan de Freitas (alandefreitas@gmail.com)

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

    Official repository: https://github.com/alandefreitas/antora-cpp-tagfiles-extension
*/

'use strict'

const {createHash} = require('node:crypto')
const expandPath = require('@antora/expand-path-helper')
const fs = require('node:fs')
const {promises: fsp} = fs
const getUserCacheDir = require('cache-directory')
const git = require('isomorphic-git')
const {globStream} = require('fast-glob')
const ospath = require('node:path')
const {posix: path} = ospath
const posixify = ospath.sep === '\\' ? (p) => p.replace(/\\/g, '/') : undefined
const {pipeline, Writable} = require('node:stream')
const forEach = (write) => new Writable({objectMode: true, write})
const yaml = require('js-yaml')
const assert = require('node:assert')
const axios = require('axios')
const {spawn} = require('node:child_process')
const {PassThrough} = require('node:stream')

const GLOB_OPTS = {ignore: ['.git'], objectMode: true, onlyFiles: false, unique: false}
const PACKAGE_NAME = 'cpp-reference-extension'
const IS_WIN = process.platform === 'win32'
const DBL_QUOTE_RX = /"/g
const QUOTE_RX = /["']/

/**
 * CppReferenceExtension is an extension that includes an extra module in the Antora
 * components with reference pages for C++ libraries and tools.
 *
 * The configuration files be registered in the components as ext.cpp-reference.config.
 * The playbook can include dependencies for the C++ library whose directories
 * can be accessed with environment variables.
 *
 * The class registers itself to the generator context and listens to the
 * 'contentAggregated' event.
 *
 * When this event is triggered, the class reads sets up the environment described in the
 * playbook and creates the reference pages for each component that defines the
 * configuration file.
 *
 * See https://docs.antora.org/antora/latest/extend/class-based-extension/
 *
 * @class
 * @property {Object} context - The generator context.
 * @property {Array} tagfiles - An array of tagfile objects.
 * @property {Object} logger - The logger object.
 * @property {Object} config - The configuration object.
 * @property {Object} playbook - The playbook object.
 *
 */
class CppReferenceExtension {
    static register({config, playbook}) {
        new CppReferenceExtension(this, {config, playbook})
    }

    constructor(generatorContext, {config, playbook}) {
        this.context = generatorContext
        const onContentAggregatedFn = this.onContentAggregated.bind(this)
        this.context.once('contentAggregated', onContentAggregatedFn)

        this.MrDocsExecutable = undefined

        // https://www.npmjs.com/package/@antora/logger
        // https://github.com/pinojs/pino/blob/main/docs/api.md
        this.logger = this.context.getLogger(PACKAGE_NAME)
        this.logger.debug('Registering cpp-reference-extension')

        this.config = config
        this.createWorktrees = config.createWorktrees || 'auto'
        // playbook = playbook
    }

    /**
     * Event handler for the 'contentAggregated' event.
     *
     * This event is triggered after all the content sources have been cloned and
     * the aggregate of all the content has been created.
     *
     * This method reads the component options and parses them.
     * The reference is generated for each component that defines the
     * reference options.
     *
     * @param {Object} playbook - The playbook object.
     * @param {Object} siteAsciiDocConfig - The AsciiDoc configuration for the site.
     * @param {Object} siteCatalog - The site catalog object.
     * @param {Array} contentAggregate - The aggregate of all the content.
     */
    async onContentAggregated({playbook, siteAsciiDocConfig, siteCatalog, contentAggregate}) {
        this.logger.debug('Reading component options')
        this.logger.debug(CppReferenceExtension.objectSummary(this.config), 'Config')
        this.logger.debug(CppReferenceExtension.objectSummary(playbook), 'Playbook')
        this.findCXXCompilers()
        const {cacheDir, gitCache, managedWorktrees} = await this.initializeWorktreeManagement(playbook);
        await this.setupDependencies(playbook, cacheDir)
        await this.setupMrDocs(playbook, cacheDir)
        for (const componentVersionBucket of contentAggregate.slice()) {
            await this.processComponentVersionBucket(componentVersionBucket, playbook, cacheDir, gitCache, managedWorktrees)
        }
        await this.performWorktreeRemovals(managedWorktrees)
    }

    /**
     * Finds and sets the environment variables for C++ and C compilers.
     *
     * The method first tries to find each compiler executable in the system's PATH. If it doesn't find it,
     * it then checks the input environment variables for a path. If a compiler is still not found, it
     * throws an error and exits the process.
     *
     * If a compiler is found, the method sets the output environment variables with the compiler path.
     *
     * @throws {Error} If a compiler is not found in the system's PATH or the input environment variables.
     */
    findCXXCompilers() {
        /*
            Find C++ compilers
            Clang++ is preferred over g++ and cl for MrDocs
        */
        const compilerConfigs = [
            {
                title: 'C++',
                executableNames: ['clang++', 'g++', 'cl'],
                inputEnv: ['CXX_COMPILER', 'CXX'],
                outputEnv: ['CMAKE_CXX_COMPILER', 'CXX']
            },
            {
                title: 'C',
                executableNames: ['clang', 'gcc', 'cl'],
                inputEnv: ['C_COMPILER', 'CC'],
                outputEnv: ['CMAKE_C_COMPILER', 'CC']
            }
        ]

        for (const {title, executableNames, inputEnv, outputEnv} of compilerConfigs) {
            let compilerExecutable = CppReferenceExtension.findExecutable(executableNames)
            if (!compilerExecutable) {
                for (const envVar of inputEnv) {
                    compilerExecutable = process.env[envVar]
                    if (compilerExecutable) {
                        break
                    }
                }
            }
            if (compilerExecutable && process.platform === "win32") {
                compilerExecutable = compilerExecutable.replace(/\\/g, '/')
            }
            if (compilerExecutable === undefined) {
                this.logger.error(`Could not find a ${title} compiler. Please set the ${inputEnv[0]} environment variable.`)
                process.exit(1)
            }
            for (const envVar of outputEnv) {
                process.env[envVar] = compilerExecutable
            }
            const compilerBasename = path.basename(compilerExecutable).replace(/\.exe$/, '')
            this.logger.debug(`${title} compiler: ${compilerBasename} (${compilerExecutable})`)
        }
    }

    /**
     * Sets up the dependencies for the playbook.
     *
     * This method iterates over the dependencies defined in the configuration. Each dependency
     * includes the name, repository URL, optional tag, and environment variable.
     *
     * The method first checks if the required fields ('name', 'repo', 'variable') are present in the dependency.
     * If a required field is missing, it logs an error and skips the dependency.
     *
     * If all required fields are present, the method checks if the dependency already exists in
     * the dependencies directory.
     * If it does, it logs a message and skips the cloning process.
     *
     * If it doesn't, it clones the repository to the dependencies directory.
     *
     * If a tag is specified, it clones the specific tag.
     *
     * After cloning, it updates the submodules.
     *
     * Finally, the method sets the environment variable with the path to the cloned repository.
     *
     * @param {Object} playbook - The playbook object.
     * @param {string} cacheDir - The cache directory.
     * @throws {Error} If a required field is missing in a dependency.
     */
    async setupDependencies(playbook, cacheDir) {
        const DependenciesDir = path.join(cacheDir, 'dependencies')
        const dependencies = this.config.dependencies
        for (const dependency of dependencies) {
            // Check if dependency has a name
            if (dependency.name === undefined) {
                if (typeof dependency.repo === 'string') {
                    this.logger.error(`Dependency name is required for ${dependency.repo}`)
                } else {
                    this.logger.error('Dependency name is required')
                }
                continue
            }

            const requiredFields = ['name', 'repo', 'variable']
            let missingRequired = false
            for (const field of requiredFields) {
                if (!dependency[field]) {
                    this.logger.error(`Dependency field "${field}" is required for ${dependency.name}`)
                    missingRequired = true
                }
            }
            if (missingRequired) {
                continue
            }

            const {name, repo, tag, variable, systemEnv, cloneSubmodules} = dependency
            if (!name) {
                this.logger.error(`Dependency name is required (${repo})`)
                continue
            }
            if (!repo) {
                this.logger.error(`Dependency repo is required (${name})`)
                continue
            }
            if (!variable) {
                this.logger.error(`Dependency variable is required (${name})`)
                continue
            }

            // Check if the dependency is already available from systemEnv
            if (systemEnv && process.env[variable]) {
                // Check if this is a directory that exists
                const dependencyPath = process.env[variable]
                const dependencyPathExists = await CppReferenceExtension.fileExists(dependencyPath)
                const dependencyPathIsDir = dependencyPathExists && (await fsp.stat(dependencyPath)).isDirectory()
                if (dependencyPathIsDir) {
                    this.logger.debug(`Dependency ${name} already exists at ${dependencyPath}`)
                    if (variable !== systemEnv) {
                        process.env[variable] = dependencyPath
                    }
                    continue
                }
            }

            let cloneDir = path.join(DependenciesDir, name)
            if (tag) {
                cloneDir = path.join(cloneDir, tag)
            }
            if (await CppReferenceExtension.fileExists(cloneDir)) {
                this.logger.debug(`Dependency ${name} already exists at ${cloneDir}`)
            } else {
                this.logger.debug(`Cloning ${repo} to ${cloneDir}`)
                await this.runCommand('git', ['clone', repo, '--depth', '1', ...(tag ? ['--branch', tag] : []), cloneDir])
                const skipCloningSubmodules = cloneSubmodules === false || cloneSubmodules === 'false'
                if (!skipCloningSubmodules) {
                    await this.runCommand('git', ['submodule', 'update', '--init', '--recursive'], {cwd: cloneDir})
                }
            }
            process.env[variable] = cloneDir
        }
    }

    /**
     * Sets up MrDocs for the playbook.
     *
     * This method downloads the latest release of MrDocs from the GitHub repository, extracts it,
     * and sets the environment variables.
     *
     * MrDocs is a tool used to generate C++ reference documentation.
     *
     * The method first determines the directories for the playbook, build, generated files, and MrDocs tree.
     * It then sends a GET request to the GitHub API to get the releases of MrDocs.
     *
     * The method iterates over the releases and finds the download URL for the latest release
     * that has binaries for the current platform.
     * If no such release is found, it logs an error and exits the process.
     *
     * If a release is found, the method determines the download directory, release tag name,
     * version subdirectory, and MrDocs executable path.
     *
     * If the MrDocs executable already exists, it logs a message and skips the download process.
     *
     * If the MrDocs executable doesn't exist, the method downloads the release from the found URL,
     * extracts it, and removes the downloaded file.
     *
     * It then checks if the MrDocs executable exists. If it does, it sets the environment variables and the MrDocs executable path.
     * If it doesn't, it logs an error and exits the process.
     *
     * @param {Object} playbook - The playbook object.
     * @param {string} cacheDir - The cache directory.
     * @throws {Error} If the MrDocs executable is not found.
     */
    async setupMrDocs(playbook, cacheDir) {
        const mrDocsTreeDir = path.join(cacheDir, 'mrdocs')
        const releasesResponse = await axios.get('https://api.github.com/repos/cppalliance/mrdocs/releases')
        const releasesInfo = releasesResponse.data
        this.logger.debug(`Found ${releasesInfo.length} MrDocs releases`)
        let downloadUrl = undefined
        let downloadRelease = undefined
        for (const latestRelease of releasesInfo) {
            this.logger.debug(`Latest release: ${latestRelease['tag_name']}`)
            const latestAssets = latestRelease['assets'].map(asset => asset['browser_download_url'])
            this.logger.debug(`Latest assets: ${latestAssets.join(', ')}`)
            const releaseFileSuffix = process.platform === "win32" ? 'win64.7z' : 'Linux.tar.gz'
            downloadUrl = latestAssets.find(asset => asset.endsWith(releaseFileSuffix))
            downloadRelease = latestRelease
            if (downloadUrl) {
                break
            }
            this.logger.warn(`Could not find MrDocs binaries in ${latestRelease['tag_name']} release for ${process.platform}`)
        }
        if (!downloadUrl) {
            this.logger.error(`Could not find MrDocs binaries for ${process.platform}`)
            process.exit(1)
        }
        const mrdocsDownloadDir = path.join(mrDocsTreeDir, process.platform)
        const releaseTagname = downloadRelease['tag_name']
        const versionSubdir = releaseTagname.endsWith('-release') ? releaseTagname.slice(0, -8) : downloadRelease['tag_name']
        const mrdocsExtractDir = path.join(mrdocsDownloadDir, versionSubdir)
        const platformExtension = process.platform === 'win32' ? '.exe' : ''
        const mrdocsExecPath = path.join(mrdocsExtractDir, 'bin', 'mrdocs') + platformExtension

        if (await CppReferenceExtension.fileExists(mrdocsExecPath)) {
            this.logger.debug(`MrDocs already exists at ${mrdocsExtractDir}`)
        } else {
            const downloadFilename = path.basename(downloadUrl)
            const downloadPath = path.join(mrdocsDownloadDir, downloadFilename)
            console.log(`Downloading ${downloadUrl} to ${mrdocsDownloadDir}...`)
            await this.downloadAndDecompress(downloadUrl, downloadPath, mrdocsExtractDir)
            await fsp.rm(downloadPath)
            console.log(`Extracted ${downloadFilename} to ${mrdocsExtractDir}`)
        }
        this.logger.debug(`Looking for MrDocs executable at ${mrdocsExecPath}`)
        if (await CppReferenceExtension.fileExists(mrdocsExecPath)) {
            this.logger.debug(`Found MrDocs executable at ${mrdocsExecPath}`)
            process.env.MRDOCS_ROOT = mrdocsExtractDir
            process.env.PATH = `${process.env.PATH}${path.delimiter}${path.join(mrdocsExtractDir, 'bin')}`
            this.MrDocsExecutable = mrdocsExecPath
        } else {
            this.logger.error(`Could not find MrDocs executable at ${mrdocsExecPath}`)
            process.exit(1)
        }
    }

    /**
     * Initializes the worktree management
     *
     * A worktree in Git is a separate working copy of the same repository
     * allowing you to work on two different branches at the same time.
     *
     * This method is used to initialize the worktree management by determining
     * the cache directory and creating it if it doesn't exist,
     * initializing a cache for git repositories, and initializing
     * a map to manage worktrees.
     *
     * It takes one argument: `playbook`. `playbook` is the playbook object.
     *
     * The method returns a promise that resolves with an object
     * containing `cacheDir`, `gitCache`, and `managedWorktrees`.
     * `cacheDir` is the determined cache directory. `gitCache`
     * is the initialized cache for git repositories. `managedWorktrees`
     * is the initialized map to manage worktrees.
     *
     * @param {Object} playbook - The playbook object.
     * @returns {Promise<Object>} A promise that resolves with an
     * object containing `cacheDir`, `gitCache`, and `managedWorktrees`.
     */
    async initializeWorktreeManagement(playbook) {
        // Determine the cache directory and create it if it doesn't exist
        const cacheDir = ospath.join(CppReferenceExtension.getBaseCacheDir(playbook), 'reference-collector')
        await fsp.mkdir(cacheDir, {recursive: true})
        this.logger.debug(`Cache directory: ${cacheDir}`)

        // Initialize a cache for git repositories
        const gitCache = {}

        // Initialize a map to manage worktrees
        const managedWorktrees = new Map()

        return {cacheDir, gitCache, managedWorktrees};
    }

    /**
     * Processes a component version bucket in the content aggregate.
     *
     * A component version bucket is an object that contains the component name,
     * version, title, startPage, asciidoc, nav, ext, files, and origins.
     * "origins" contains various versions of the component.
     *
     * This method is used to process a component version bucket in the
     * content aggregate.
     *
     * @param {Object} componentVersionBucket - The component version bucket to be processed.
     * @param {Object} playbook - The playbook object.
     * @param {string} cacheDir - The cache directory.
     * @param {Object} gitCache - The cache for git repositories.
     * @param {Map} managedWorktrees - The map to manage worktrees.
     */
    async processComponentVersionBucket(componentVersionBucket, playbook, cacheDir, gitCache, managedWorktrees) {
        this.logger.debug(CppReferenceExtension.objectSummary(componentVersionBucket), 'Process component version bucket')
        for (const origin of componentVersionBucket.origins) {
            await this.processOrigin(origin, componentVersionBucket, playbook, cacheDir, gitCache, managedWorktrees);
        }
    }

    /**
     * Processes an origin of a component version bucket.
     *
     * An origin is an object that contains the type (e.g. git), url,
     * gitdir (e.g. path/to/repo/.git), reftype (e.g. branch),
     * refname (e.g. develop), branch (e.g. develop), startPath (e.g. docs),
     * worktree, fileUriPattern, webUrl, editUrlPattern, and descriptor.
     *
     * The descriptor is an object that contains the contents of the antora.yml file
     * such as name, version, title, startPage, asciidoc, nav, and ext.
     *
     * This method is used to process an origin of a component version bucket.
     * It determines the directory for the worktree, creates a context for
     * expanding paths, and creates a list of normalized collectors from the collector
     * configuration.
     *
     * @param {Object} origin - The origin to be processed.
     * @param {Object} componentVersionBucket - The component version bucket that contains the origin.
     * @param {Object} playbook - The playbook object.
     * @param {string} cacheDir - The cache directory.
     * @param {Object} gitCache - The cache for git repositories.
     * @param {Map} managedWorktrees - The map to managed worktrees.
     */
    async processOrigin(origin, componentVersionBucket, playbook, cacheDir, gitCache, managedWorktrees) {
        const {name, version} = componentVersionBucket
        const {url, gitdir, refname, reftype, remote, worktree, startPath, descriptor} = origin
        this.logger.debug(`Processing origin ${url || gitdir} (${reftype}: ${refname}) at path ${startPath}`)
        this.logger.debug(CppReferenceExtension.objectSummary(origin), 'Origin')

        // Get the reference collector configuration from the descriptor
        // The reference collector configuration is an array of objects
        // because components are allowed to define multiple collectors
        let collectorConfigs = descriptor?.ext?.cppReference || []
        if (!Array.isArray(collectorConfigs)) {
            collectorConfigs = [collectorConfigs]
        }
        if (!collectorConfigs.length) {
            this.logger.warn(`No reference collector configuration found for component ${name} version ${version}`)
            return
        }
        this.logger.debug(CppReferenceExtension.objectSummary(collectorConfigs), 'Component C++ reference configuration')

        // Determine the directory for the worktree.
        // A worktree is a Git feature that allows you to have multiple
        // branches of a repository checked out at once.
        // Each worktree has its own directory.
        let worktreeDir = worktree
        let worktreeConfig = collectorConfigs[0].worktree || {}
        if (!worktreeConfig) {
            worktreeConfig = worktreeConfig === false ? {create: 'always'} : {}
        }
        const createWorktree =
            !worktree || ('create' in worktreeConfig ? worktreeConfig.create : this.createWorktrees) === 'always'
        const checkoutWorktree = worktreeConfig.checkout !== false
        if (createWorktree) {
            this.logger.debug(`Worktree directory not provided for ${name} version ${version}`)
            worktreeDir = await this.setupManagedWorktree(worktreeConfig, checkoutWorktree, origin, cacheDir, managedWorktrees);
        }

        // Store the worktree directory in the origin
        origin.collectorWorktree = worktreeDir
        this.logger.debug(`Worktree directory: ${worktreeDir}`)

        // Create a context for expanding paths
        const expandPathContext = {
            // The base directory for expanding paths
            base: worktreeDir,
            // The current working directory
            cwd: worktreeDir,
            // The start path for expanding paths
            dot: ospath.join(worktreeDir, startPath)
        }
        this.logger.debug(expandPathContext, 'Expand path context')

        // If the worktree doesn't exist, either checkout the worktree or create the directory
        if (createWorktree) {
            if (checkoutWorktree) {
                this.logger.debug(`Checking out worktree: ${worktreeDir}`)
                const cache = gitCache[gitdir] || (gitCache[gitdir] = {})
                const ref = `refs/${reftype === 'branch' ? 'head' : reftype}s/${refname}`
                this.logger.debug(cache, 'Cache')
                this.logger.debug(`Ref: ${ref}.`)
                await this.prepareWorktree({
                    fs,
                    cache,
                    dir: worktreeDir,
                    gitdir,
                    ref,
                    remote,
                    bare: worktree === undefined
                })
                this.logger.debug(`Checked out worktree: ${worktreeDir}`)
            } else {
                this.logger.debug(`Creating worktree directory: ${worktreeDir}`)
                await fsp.mkdir(worktreeDir, {recursive: true})
            }
        } else {
            this.logger.debug(`Using existing worktree directory: ${worktreeDir}`)
        }

        // Create a list of normalized collectors from the collector configuration
        let collectors = []
        for (const collectorConfig of collectorConfigs) {
            // If a config file is specified in the collectorConfig configuration, check if it exists
            let mrdocsConfigFile
            const defaultMrDocsConfigLocations = ['mrdocs.yml', 'docs/mrdocs.yml', 'doc/mrdocs.yml'];
            const mrdocsConfigCandidates = [collectorConfig.config].concat(defaultMrDocsConfigLocations)
            this.logger.debug(`Looking for mrdocs.yml file in ${startPath} at locations ${mrdocsConfigCandidates.join(', ')} for component ${name} version ${version}`)
            for (const candidate of mrdocsConfigCandidates) {
                if (candidate) {
                    const candidateBasePaths = [expandPathContext.base, expandPathContext.dot]
                    this.logger.debug(`Base paths: ${candidateBasePaths.join(', ')}`)
                    for (const basePath of candidateBasePaths) {
                        const candidatePath = path.join(basePath, candidate)
                        this.logger.debug(`Checking candidate path: ${candidatePath}`)
                        if (fs.existsSync(candidatePath)) {
                            mrdocsConfigFile = candidatePath
                            this.logger.debug(`Found mrdocs.yml file: ${mrdocsConfigFile}`)
                            break
                        }
                    }
                }
                if (mrdocsConfigFile) {
                    break
                }
            }
            if (!mrdocsConfigFile) {
                this.logger.warn(`No mrdocs.yml file found in ${startPath} at locations ${mrdocsConfigCandidates.join(', ')} for component ${name} version ${version}`)
                continue
            }
            this.logger.debug(`Using mrdocs.yml file: ${mrdocsConfigFile}`)

            collectors.push({
                config: mrdocsConfigFile
            })
        }
        this.logger.debug(collectors, 'Collectors')

        // For each collector, perform clean, run, and scan operations
        for (const collector of collectors) {
            await this.processCollector(collector, origin, componentVersionBucket, playbook, worktreeDir, worktree, cacheDir);
        }
    }

    /**
     * Prepares a worktree directory for a given origin.
     *
     * A worktree in Git is a separate working copy of the same repository
     * allowing you to work on two different branches at the same time.
     *
     * This method is used to prepare a worktree directory for a given origin.
     * It determines whether to check out the worktree and whether to keep the
     * worktree after use.
     *
     * It generates a name for the worktree directory and checks if the
     * worktree directory is already being managed.
     *
     * - If the worktree directory is already being managed, it adds the origin to it.
     * - If the worktree directory is not being managed, it adds it to the managed worktrees map.
     * - If the worktree is not being checked out or it's not being kept, it removes the directory.
     *
     * @param {Object} worktreeConfig - The worktree configuration object.
     * @param {boolean} checkoutWorktree - Whether to checkout the worktree.
     * @param {Object} origin - The origin object.
     * @param {string} cacheDir - The cache directory.
     * @param {Map} managedWorktrees - The map to manage worktrees.
     * @returns {Promise<string>} A promise that resolves with the worktree directory.
     */
    async setupManagedWorktree(worktreeConfig, checkoutWorktree, origin, cacheDir, managedWorktrees) {
        // Determine whether we should keep the worktree after use.
        // By default, we don't keep it unless explicitly set to true.
        const keepWorktree =
            'keep' in worktreeConfig ?
                worktreeConfig.keep :
                'keepWorktrees' in this.config ?
                    this.config.keepWorktrees === true :
                    false
        this.logger.debug(`Creating worktree for ${origin.url} with keepWorktree=${keepWorktree}`)

        // Generate a name for the worktree directory and join it with
        // the cache directory path.
        const worktreeFolderName = CppReferenceExtension.generateWorktreeFolderName(origin, keepWorktree);
        let worktreeDir = ospath.join(cacheDir, 'worktrees', worktreeFolderName)
        this.logger.debug(`Worktree directory: ${worktreeDir}`)

        // Check if the worktree directory is already being managed.
        // If it is, we add the origin to it.
        // Otherwise, we create a new entry in the managed worktrees map.
        if (managedWorktrees.has(worktreeDir)) {
            this.logger.debug(`Worktree directory ${worktreeDir} is already being managed`)
            managedWorktrees.get(worktreeDir).origins.add(origin)
            // If we're not checking out the worktree, we remove the directory.
            if (!checkoutWorktree) {
                this.logger.debug(`Removing worktree directory ${worktreeDir} as we're not checking it out`)
                await fsp.rm(worktreeDir, {force: true, recursive: true})
            }
        } else {
            this.logger.debug(`Worktree directory ${worktreeDir} is not being managed`)
            // If the worktree directory is not being managed, we add it to the managed worktrees map.
            managedWorktrees.set(worktreeDir, {origins: new Set([origin]), keep: keepWorktree})
            // If we're not checking out the worktree, or we're not keeping it, we remove the directory.
            if (!checkoutWorktree || keepWorktree !== true) {
                this.logger.debug(`Removing worktree directory ${worktreeDir} as we're not checking it out or keeping it`)
                await fsp.rm(worktreeDir, {
                    force: true,
                    recursive: true
                })
            }
        }
        return worktreeDir;
    }

    /**
     * Prepares a git worktree from the specified gitdir, making use of the existing clone.
     *
     * If the worktree already exists from a previous iteration, the worktree is reset.
     *
     * A valid worktree is one that contains a .git/index file.
     * Otherwise, a fresh worktree is created.
     *
     * If the gitdir contains an index file, that index file is temporarily overwritten to
     * prepare the worktree and later restored before the function returns.
     *
     * @param {Object} repo - The repository object.
     */
    async prepareWorktree(repo) {
        const {dir: worktreeDir, gitdir, ref, remote = 'origin', bare, cache} = repo
        this.logger.debug(`Preparing worktree for ${worktreeDir} from ${gitdir} at ${ref}`)
        delete repo.remote
        const currentIndexPath = ospath.join(gitdir, 'index')
        this.logger.debug(`Current index: ${currentIndexPath}`)
        const currentIndexPathBak = currentIndexPath + '~'
        this.logger.debug(`Current index backup: ${currentIndexPathBak}`)
        const restoreIndex = (await fsp.rename(currentIndexPath, currentIndexPathBak).catch(() => false)) === undefined
        this.logger.debug(`Restore index: ${restoreIndex} because it was not possible to rename ${currentIndexPath} to ${currentIndexPathBak}`)
        const worktreeGitdir = ospath.join(worktreeDir, '.git')
        this.logger.debug(`Worktree gitdir: ${worktreeGitdir}`)
        const worktreeIndexPath = ospath.join(worktreeGitdir, 'index')
        this.logger.debug(`Worktree index: ${worktreeIndexPath}`)
        try {
            let force = true
            try {
                await CppReferenceExtension.mv(worktreeIndexPath, currentIndexPath)
                this.logger.debug(`Moved ${worktreeIndexPath} to ${currentIndexPath}`)
                await CppReferenceExtension.removeUntrackedFiles(repo)
                this.logger.debug(`Removed untracked files from ${worktreeDir}`)
            } catch {
                this.logger.debug(`Could not move ${worktreeIndexPath} to ${currentIndexPath}`)
                force = false
                // index file not needed in this case
                await fsp.unlink(currentIndexPath).catch(() => undefined)
                await fsp.rm(worktreeDir, {recursive: true, force: true})
                await fsp.mkdir(worktreeGitdir, {recursive: true})
                this.logger.debug(`Created worktree directory ${worktreeDir}`)
                Reflect.ownKeys(cache).forEach((it) => it.toString() === 'Symbol(PackfileCache)' || delete cache[it])
                this.logger.debug(`Removed cache for ${worktreeDir}`)
            }
            let head
            if (ref.startsWith('refs/heads/')) {
                head = `ref: ${ref}`
                const branchName = ref.slice(11)
                if (bare || !(await git.listBranches(repo)).includes(branchName)) {
                    await git.branch({
                        ...repo,
                        ref: branchName,
                        object: `refs/remotes/${remote}/${branchName}`,
                        force: true
                    })
                }
            } else {
                head = await git.resolveRef(repo)
            }
            this.logger.debug(`Checking out HEAD: ${head}`)
            await git.checkout({...repo, force, noUpdateHead: true, track: false})
            this.logger.debug(`Checked out HEAD: ${head} at ${worktreeDir}`)
            await fsp.writeFile(ospath.join(worktreeGitdir, 'commondir'), `${gitdir}\n`, 'utf8')
            this.logger.debug(`Wrote commondir: ${gitdir}`)
            const headPath = ospath.join(worktreeGitdir, 'HEAD');
            await fsp.writeFile(headPath, `${head}\n`, 'utf8')
            this.logger.debug(`Wrote HEAD path: ${headPath}`)
            await CppReferenceExtension.mv(currentIndexPath, worktreeIndexPath)
            this.logger.debug(`Moved ${currentIndexPath} to ${worktreeIndexPath}`)
        } finally {
            if (restoreIndex) await fsp.rename(currentIndexPathBak, currentIndexPath)
        }
    }

    static mv(from, to) {
        return fsp.cp(from, to).then(() => fsp.rm(from))
    }

    static removeUntrackedFiles(repo) {
        const trees = [git.STAGE({}), git.WORKDIR()]
        const map = (relpath, [sEntry]) => {
            if (relpath === '.') return
            if (relpath === '.git') return null
            if (sEntry == null) return fsp.rm(ospath.join(repo.dir, relpath), {recursive: true}).then(invariably.null)
            return sEntry.mode().then((mode) => (mode === 0o120000 ? null : undefined))
        }
        return git.walk({...repo, trees, map})
    }

    /**
     * Processes a collector of a component version bucket.
     *
     * A collector is an object that contains the configuration for
     * generating the C++ reference documentation.
     *
     * This method is used to process a collector of a component
     * version bucket. It determines the directory for the reference,
     * ensures the reference output directory exists and is clean,
     * sets up MrDocs, and creates an index.adoc file in the
     * reference output directory.
     *
     * It then scans the directories for the reference output
     * and adds them to the target files.
     *
     * @param {Object} collector - The collector to be processed.
     * @param {Object} origin - The origin object.
     * @param {Object} componentVersionBucket - The component version bucket that contains the collector.
     * @param {Object} playbook - The playbook object.
     * @param {string} worktreeDir - The worktree directory.
     * @param {string} cacheDir - The cache directory.
     * @param worktree
     */
    async processCollector(collector, origin, componentVersionBucket, playbook, worktreeDir, worktree, cacheDir) {
        const {name, title, version = []} = componentVersionBucket;
        this.logger.debug(`Processing collector for ${title} (${name}) version ${version}`)

        // Determine the directory for the reference
        assert(typeof (version) === 'string', 'Version should be a string')
        const referenceOutputDir =
            (version && typeof version === 'string') ?
                path.join(cacheDir, 'reference', name, 'versioned', version) :
                path.join(cacheDir, 'reference', name, 'main')
        this.logger.debug(`Reference output directory: ${referenceOutputDir}`)

        // Make sure the reference output directory exists and it's clean
        if (fs.existsSync(referenceOutputDir)) {
            await fsp.rm(referenceOutputDir, {recursive: true, force: true})
        }
        await fsp.mkdir(referenceOutputDir, {recursive: true})

        // Generate reference documentation with MrDocs
        await this.runCommand(this.MrDocsExecutable, [
            `--config=${collector.config}`,
            `--output=${referenceOutputDir}`,
            `--generate=adoc`,
            `--multipage=true`
        ], {cwd: worktreeDir, quiet: playbook.runtime?.quiet})

        await this.scanDirectories(referenceOutputDir, origin, componentVersionBucket, worktreeDir, worktree);
    }

    /**
     * Scans directories and adds files to the target files.
     *
     * This method is used to scan directories and add files to the target files.
     *
     * It determines the module path for each file and checks if the file already exists in the target files.
     * If the file exists, it updates the contents and stats of the file in the target files.
     * If the file doesn't exist, it adds the file to the target files.
     *
     * The method does not return a value.
     *
     * @param {string} referenceOutputDir - The directory where the reference output is stored.
     * @param {Object} origin - The origin object.
     * @param {Object} componentVersionBucket - The component version bucket that contains the files.
     * @param {string} worktreeDir - The worktree directory.
     * @param {string} worktree - The original worktree directory.
     */
    async scanDirectories(referenceOutputDir, origin, componentVersionBucket, worktreeDir, worktree) {
        // Scan directories
        const referenceModuleName = 'reference'
        const relModulePrefix = path.join('modules', referenceModuleName, 'pages')
        const files = await CppReferenceExtension.srcFs(referenceOutputDir, '**/*')
        const targetFiles = componentVersionBucket.files
        for (const file of files) {
            this.logger.debug(CppReferenceExtension.objectSummary(file), 'Scanning File')
            const relpath = file.path
            const modulePath = path.join(relModulePrefix, relpath)
            const existingFile = targetFiles.find((it) => it.path === modulePath)
            if (existingFile) {
                Object.assign(existingFile, {contents: file.contents, stat: file.stat})
            } else {
                Object.assign(file, {path: path.join(relModulePrefix, file.path)})
                Object.assign(file.src, {
                    path: path.join(relModulePrefix, file.src.path),
                    abspath: path.join(modulePath, file.src.path)
                })
                this.logger.debug(CppReferenceExtension.objectSummary(file), `Adding reference file to ${modulePath}`)
                const src = file.src
                const scannedRelpath = src.abspath.slice(worktreeDir.length + 1)
                Object.assign(src, {
                    origin,
                    scanned: posixify ? posixify(scannedRelpath) : scannedRelpath
                })
                if (!worktree) {
                    Object.assign(src, {realpath: src.abspath, abspath: src.scanned})
                }
                targetFiles.push(file)
            }
        }
    }

    /**
     * Performs the removal of worktrees.
     *
     * A worktree in Git is a separate working copy of the same repository
     * allowing you to work on two different branches at the same time.
     *
     * This method is used to perform the removal of worktrees.
     * It prepares for deferred worktree removals
     * and then performs the deferred worktree removals.
     *
     * @param {Map} managedWorktrees - The map of worktrees that are being managed.
     */
    async performWorktreeRemovals(managedWorktrees) {
        const deferredWorktreeRemovals = await this.prepareDeferredWorktreeRemovals(managedWorktrees);
        await this.performDeferredWorktreeRemovals(deferredWorktreeRemovals);
    }

    /**
     * Prepares for the deferred removal of worktrees.
     *
     * A worktree in Git is a separate working copy of the same repository
     * allowing you to work on two different branches at the same time.
     *
     * This method is used to prepare for the deferred removal of worktrees.
     * It iterates over the managed worktrees and checks the 'keep' property of each worktree.
     * If 'keep' is true, the worktree is skipped.
     * If 'keep' is a string that starts with 'until:', the worktree removal is deferred until the specified event.
     * Otherwise, the worktree is removed immediately.
     *
     * The method returns a promise that resolves with a map of deferred worktree removals.
     * Each entry in the map is an array of worktrees to be removed when a specific event occurs.
     *
     * @param {Map} managedWorktrees - The map of worktrees that are being managed.
     * @returns {Promise<Map>} A promise that resolves with a map of deferred worktree removals.
     */
    async prepareDeferredWorktreeRemovals(managedWorktrees) {
        // Prepare for deferred worktree removals
        this.logger.debug('Preparing for manual worktree removals')
        const deferredWorktreeRemovals = new Map()
        for (const [worktreeDir, {origins, keep}] of managedWorktrees) {
            this.logger.debug(`Managing worktree directory: ${worktreeDir}`)
            if (keep === true) continue
            if (typeof keep === 'string' && keep.startsWith('until:')) {
                const eventName = keep === 'until:exit' ? 'contextClosed' : keep.slice(6)
                const removal = {worktreeDir, origins}
                const removals = deferredWorktreeRemovals.get(eventName)
                removals ? removals.push(removal) : deferredWorktreeRemovals.set(eventName, [removal])
                continue
            }
            await CppReferenceExtension.removeWorktree(worktreeDir, origins)
        }
        return deferredWorktreeRemovals
    }

    /**
     * Performs the deferred removal of worktrees.
     *
     * A worktree in Git is a separate working copy of the same repository
     * allowing you to work on two different branches at the same time.
     *
     * This method is used to perform the deferred removal of worktrees.
     * It iterates over the deferred worktree removals and for each removal,
     * it sets up an event listener that triggers the removal of the worktree
     * when the specified event occurs.
     *
     * It takes one argument: `deferredWorktreeRemovals`. `deferredWorktreeRemovals` is a map
     * where each entry is an array of worktrees to be removed when a specific event occurs.
     *
     * The method does not return a value.
     *
     * @param {Map} deferredWorktreeRemovals - A map of deferred worktree removals.
     */
    async performDeferredWorktreeRemovals(deferredWorktreeRemovals) {
        // Perform deferred worktree removals
        for (const [eventName, removals] of deferredWorktreeRemovals) {
            this.context.once(eventName, () =>
                Promise.all(
                    removals.map(({worktreeDir, origins}) =>
                        CppReferenceExtension.removeWorktree(worktreeDir, origins)))
            )
        }
    }

    /**
     * Generates a unique folder name for a worktree.
     *
     * A worktree in Git is a separate working copy of the same repository
     * allowing you to work on two different branches at the same time.
     *
     * This method generates a unique folder name for a worktree based on the
     * repository URL, the directory of the Git repository, the reference name
     * (branch name), and the worktree directory.
     *
     * @param {Object} options - An object containing `url`, `gitdir`, `refname`, and `worktree` properties.
     * @param {string} options.url - The URL of the Git repository.
     * @param {string} options.gitdir - The directory of the Git repository.
     * @param {string} options.refname - The reference name (branch name).
     * @param {string} options.worktree - The worktree directory.
     * @param {boolean} keepWorktrees - Flag indicating whether to keep worktrees.
     * @returns {string} The generated folder name for the worktree.
     */
    static generateWorktreeFolderName({url, gitdir, refname, worktree}, keepWorktrees) {
        // Create a qualifier for the reference name if worktrees are to be kept
        const refnameQualifier = keepWorktrees ? '@' + refname.replace(/[/]/g, '-') : undefined

        // If worktree is undefined, generate a folder name based on the gitdir
        if (worktree === undefined) {
            const folderName = ospath.basename(gitdir, '.git')
            if (!refnameQualifier) return folderName
            const lastHyphenIdx = folderName.lastIndexOf('-')
            return `${folderName.slice(0, lastHyphenIdx)}${refnameQualifier}${folderName.slice(lastHyphenIdx)}`
        }

        // Normalize the URL or gitdir
        let normalizedUrl = (url || gitdir).toLowerCase()
        if (posixify) normalizedUrl = posixify(normalizedUrl)
        normalizedUrl = normalizedUrl.replace(/(?:[/]?\.git|[/])$/, '')

        // Create a slug based on the normalized URL and the refname qualifier
        const slug = ospath.basename(normalizedUrl) + (refnameQualifier || '')

        // Create a hash of the normalized URL
        const hash = createHash('sha1').update(normalizedUrl).digest('hex')

        // Return the slug and hash as the folder name
        return `${slug}-${hash}`
    }

    /**
     * Determines the base directory for caching.
     *
     * This static method of the `CppReferenceExtension` class is used to determine the base directory for caching.
     * It takes an object as an argument, which has two properties: `dir` and `runtime`. `dir` is aliased as `dot`,
     * and `runtime` is an object that contains a `cacheDir` property.
     *
     * If `cacheDir` is `null`, the function tries to get the user cache directory by calling the `getUserCacheDir`
     * function with a string argument. This string argument is either 'antora' or 'antora-test', depending on whether
     * the `NODE_ENV` environment variable is set to 'test'. If `getUserCacheDir` returns a falsy value, it falls back
     * to a default directory path, which is the `.cache/antora` directory inside the `dir` directory.
     *
     * If `cacheDir` is not `null`, the function calls `expandPath` with `cacheDir` and an object containing `dir` as arguments.
     *
     * @param {Object} options - An object containing `dir` and `runtime` properties.
     * @param {string} options.dir - The directory to use as a base for caching.
     * @param {Object} options.runtime - An object containing a `cacheDir` property.
     * @param {string} options.runtime.cacheDir - The preferred cache directory.
     * @returns {string} The determined cache directory.
     */
    static getBaseCacheDir({dir: dot, runtime: {cacheDir: preferredDir}}) {
        return preferredDir == null
            ? getUserCacheDir(`antora${process.env.NODE_ENV === 'test' ? '-test' : ''}`) || ospath.join(dot, '.cache/antora')
            : expandPath(preferredDir, {dot})
    }

    /**
     * Reads files from a directory and its subdirectories based on provided glob patterns.
     *
     * This static method of the `CppReferenceExtension` class is used to read files from a directory and its subdirectories.
     * It takes three arguments: `cwd`, `globs`, and `into`. `cwd` is the current working directory. `globs` is a pattern or an array of patterns
     * matching the files to be read. `into` is a directory path where the read files will be placed.
     *
     * The method returns a promise that resolves with an array of file objects. Each file object contains the file path, contents, stats, and source information.
     *
     * @param {string} cwd - The current working directory.
     * @param {string|string[]} [globs] - The glob pattern or an array of glob patterns matching the files to be read.
     * @param {string} [into] - The directory path where the read files will be placed.
     * @returns {Promise<Array>} A promise that resolves with an array of file objects.
     */
    static srcFs(cwd, globs = '**/*', into) {
        return new Promise((resolve, reject, accum = []) =>
            // Create a pipeline with a glob stream and a writable stream
            pipeline(
                // Create a glob stream with the provided glob patterns
                globStream(globs, Object.assign({cwd}, GLOB_OPTS)),
                // Create a writable stream that processes each file matched by the glob patterns
                forEach(({path: relpath, dirent}, _, done) => {
                    // If the matched file is a directory, skip it
                    if (dirent.isDirectory()) return done()
                    // Normalize the file path
                    let relpathPosix = relpath
                    // Determine the absolute path of the file
                    const abspath = posixify ? ospath.join(cwd, (relpath = ospath.normalize(relpath))) : cwd + '/' + relpath
                    // Get the file stats
                    fsp.stat(abspath).then((stat) => {
                        // Read the file contents
                        fsp.readFile(abspath).then((contents) => {
                            // Determine the basename, extension, and stem of the file
                            const basename = ospath.basename(relpathPosix)
                            const extname = ospath.extname(relpathPosix)
                            const stem = basename.slice(0, basename.length - extname.length)
                            // If an `into` directory is provided, update the file path
                            if (into) relpathPosix = path.join('.', into, relpathPosix)
                            // Add the file to the accumulator
                            accum.push({
                                path: relpathPosix,
                                contents,
                                stat,
                                src: {path: relpathPosix, basename, stem, extname, abspath},
                            })
                            // Proceed to the next file
                            done()
                        }, done)
                    }, done)
                }),
                // Resolve the promise with the accumulated files or reject it with an error
                (err) => (err ? reject(err) : resolve(accum))
            )
        )
    }

    /**
     * Asynchronously removes a worktree directory and deletes the 'collectorWorktree' property from each origin.
     *
     * A worktree in Git is a separate working copy of the same repository
     * allowing you to work on two different branches at the same time.
     *
     * This static method of the `CppReferenceExtension` class is used to remove a worktree
     * directory and delete the 'collectorWorktree' property from each origin.
     *
     * It takes two arguments: `worktreeDir` and `origins`. `worktreeDir` is the directory
     * of the worktree to be removed. `origins` is an array of origin objects.
     *
     * The method returns a promise that resolves when the worktree directory has been
     * removed and the 'collectorWorktree' property has been deleted from each origin.
     *
     * @param {string} worktreeDir - The directory of the worktree to be removed.
     * @param {Array} origins - An array of origin objects.
     * @returns {Promise} A promise that resolves when the worktree directory has been removed and the 'collectorWorktree' property has been deleted from each origin.
     */
    static async removeWorktree(worktreeDir, origins) {
        // Iterate over each origin
        for (const origin of origins) {
            // Delete the 'collectorWorktree' property from the origin
            delete origin.collectorWorktree
        }
        // Remove the worktree directory
        await fsp.rm(worktreeDir, {recursive: true})
    }


    /**
     * Finds an executable in the system's PATH.
     *
     * This static method of the `CppReferenceExtension` class is used to find an executable in the system's PATH.
     * It takes one argument: `executableName`.
     *
     * `executableName` is the name of the executable to find or an array of possible executable names.
     * On Windows, the method checks for executables with the extensions `.exe`, `.bat`, and `.cmd`.
     *
     * It can be a string or an array of strings. If it's an array, the method returns the path
     * of the first executable found.
     *
     * The method checks if the executable exists in each directory of the PATH.
     * If the executable is not found, it searches for versioned executables.
     *
     * Versioned executables are executables with a version number in the name.
     * For instance, looking for `clang` might find `clang-12` if it's the highest version
     * available.
     *
     * The method returns the path of the executable if found, or `undefined` if not found.
     *
     * @param {string|string[]} executableName - The name of the executable to find.
     * @returns {string|undefined} The path of the executable if found, or `undefined` if not found.
     */
    static findExecutable(executableName) {
        if (Array.isArray(executableName)) {
            for (const name of executableName) {
                const result = CppReferenceExtension.findExecutable(name)
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


    /**
     * Downloads a file from a given URL and decompresses it.
     *
     * This static method of the `CppReferenceExtension` class is used to download a file from a given URL and decompress it.
     * It takes three arguments: `downloadUrl`, `downloadPath`, and `extractPath`.
     * `downloadUrl` is the URL of the file to be downloaded.
     * `downloadPath` is the path where the downloaded file will be saved.
     * `extractPath` is the path where the downloaded file will be decompressed.
     *
     * The method sends a GET request to the `downloadUrl` and writes the response body to the `downloadPath`.
     * If the downloaded file is a .7z file, it decompresses it using the 7z command.
     * If the downloaded file is a .tar.gz file, it decompresses it using the tar command.
     *
     * The method does not return a value.
     *
     * @param {string} downloadUrl - The URL of the file to be downloaded.
     * @param {string} downloadPath - The path where the downloaded file will be saved.
     * @param {string} extractPath - The path where the downloaded file will be decompressed.
     */
    async downloadAndDecompress(downloadUrl, downloadPath, extractPath) {
        // ========================================================================
        // Download
        // ========================================================================
        this.logger.debug(`Downloading ${downloadUrl} to ${downloadPath}...`);
        const response = await axios.get(downloadUrl, {responseType: 'arraybuffer'});
        if (response.status !== 200) {
            this.logger.error(`Failed to download ${downloadUrl} (Error ${response.status})`)
            process.exit(1)
        }
        this.logger.debug(`Downloaded ${downloadUrl}.`);
        const fileData = Buffer.from(response.data, 'binary');
        await fsp.mkdir(path.dirname(downloadPath), {recursive: true})
        await fsp.writeFile(downloadPath, fileData);
        this.logger.debug(`File ${downloadUrl} downloaded successfully to ${downloadPath}.`);

        // ========================================================================
        // Deflate
        // ========================================================================
        await fsp.mkdir(extractPath, {recursive: true})
        const tempExtractPath = extractPath + '-temp'
        if (await CppReferenceExtension.fileExists(tempExtractPath)) {
            await fsp.rm(tempExtractPath, {recursive: true})
        }
        await fsp.mkdir(tempExtractPath, {recursive: true})
        if (path.extname(downloadPath) === '.7z') {
            await this.runCommand('7z', ['x', downloadPath, `-o${tempExtractPath}`], {output: true})
        } else if (/\.tar\.gz$/.test(downloadPath)) {
            await this.runCommand('tar', ['-vxzf', downloadPath, '-C', tempExtractPath], {output: true})
        }
        const files = await fsp.readdir(tempExtractPath);
        const nFiles = files.length
        if (nFiles === 1 && (await fsp.stat(path.join(tempExtractPath, files[0]))).isDirectory()) {
            await fsp.rm(extractPath, {recursive: true})
            await fsp.rename(path.join(tempExtractPath, files[0]), extractPath)
            await fsp.rm(tempExtractPath, {recursive: true})
        } else {
            await fsp.rm(extractPath, {recursive: true})
            await fsp.rename(tempExtractPath, extractPath)
        }
        this.logger.debug(`File decompressed successfully to ${extractPath}.`);
    }

    /**
     * Executes a command in a child process.
     *
     * This method is used to execute a command in a child process.
     *
     * The method returns a promise that resolves with the standard output of the
     * command if the command is executed successfully.
     *
     * If the command execution fails, the promise is rejected with an error.
     *
     * @param {string} cmd - The command to be executed.
     * @param {Array} [argv=[]] - The arguments to be passed to the command.
     * @param {Object} [opts={}] - The options for the command execution.
     * @param {Buffer|string} [opts.input] - The input to be passed to the command.
     * @param {boolean} [opts.output] - Flag indicating whether to output the command's standard output.
     * @param {boolean} [opts.quiet] - Flag indicating whether to suppress the command's standard output.
     * @param {boolean} [opts.implicitStdin] - Flag indicating whether to implicitly pass the standard input to the command.
     * @param {boolean} [opts.local] - Flag indicating whether to execute the command in the local directory.
     * @returns {Promise<Buffer|string>} A promise that resolves with the standard output of the command.
     * @throws {Error} If the command execution fails.
     */
    async runCommand(cmd, argv = [], opts = {}) {
        this.logger.debug({cmd, argv, opts}, 'Running command')
        if (!cmd) {
            throw new TypeError('Command not specified')
        }
        let cmdv = CppReferenceExtension.parseCommand(cmd)
        const {input, output, quiet, implicitStdin, local, ...spawnOpts} = opts
        if (input) {
            input instanceof Buffer ? implicitStdin || argv.push('-') : argv.push(input)
        }
        if (IS_WIN) {
            if (local && !cmdv[0].endsWith('.bat')) {
                const cmd0 = `${cmdv[0]}.bat`
                const cmdExists = await CppReferenceExtension.fileExists(ospath.join(opts.cwd || '', cmd0));
                if (cmdExists) {
                    cmdv[0] = cmd0
                }
            }
            cmdv = cmdv.map(CppReferenceExtension.winShellEscape)
            argv = argv.map(CppReferenceExtension.winShellEscape)
            Object.assign(spawnOpts, {shell: true, windowsHide: true})
        } else if (local) {
            cmdv[0] = `./${cmdv[0]}`
        }
        return new Promise((resolve, reject) => {
            const stdout = []
            const stderr = []
            const ps = spawn(cmdv[0], [...cmdv.slice(1), ...argv], spawnOpts)
            ps.on('close', (code) => {
                if (code === 0) {
                    if (stderr.length) {
                        process.stderr.write(stderr.join(''))
                    }
                    if (output) {
                        // adapted from https://github.com/jpommerening/node-lazystream/blob/master/lib/lazystream.js | license: MIT
                        class LazyReadable extends PassThrough {
                            constructor(fn, options) {
                                super(options)
                                this._read = function () {
                                    delete this._read // restores original method
                                    fn.call(this, options).on('error', this.emit.bind(this, 'error')).pipe(this)
                                    return this._read.apply(this, arguments)
                                }
                                this.emit('readable')
                            }
                        }

                        output === true ? resolve() : resolve(new LazyReadable(() => fs.createReadStream(output)))
                    } else {
                        resolve(Buffer.from(stdout.join('')))
                    }
                } else {
                    let msg = `Command failed with exit code ${code}: ${ps.spawnargs.join(' ')}`
                    if (stderr.length) msg += '\n' + stderr.join('')
                    reject(new Error(msg))
                }
            })
            ps.on('error', (err) => reject(err.code === 'ENOENT' ? new Error(`Command not found: ${cmdv.join(' ')}`) : err))
            ps.stdout.on('data', (data) => (output ? !quiet && process.stdout.write(data) : stdout.push(data)))
            ps.stderr.on('data', (data) => stderr.push(data))
            try {
                input instanceof Buffer ? ps.stdin.end(input) : ps.stdin.end()
            } catch (err) {
                reject(err)
            } finally {
                ps.stdin.end()
            }
        })
    }

    /**
     * Checks if a file or directory exists at the given path.
     *
     * This method uses the fs.promises API's access method to check
     * if a file or directory exists.
     *
     * If the file or directory exists, the method resolves the
     * promise with true.
     *
     * If the file or directory does not exist, the method
     * resolves the promise with false.
     *
     * @param {string} p - The path to the file or directory.
     * @returns {Promise<boolean>} A promise that resolves with true if the file
     * or directory exists, or false if it does not exist.
     */
    static fileExists(p) {
        return fsp.access(p).then(
            () => true,
            () => false
        )
    }

    /**
     * Escapes a string for use in a Windows shell command.
     *
     * This method is used to escape a string for use in a Windows shell command.
     *
     * The method checks if the first character of the string is a hyphen (-).
     * If it is, the method returns `val` as is.
     *
     * If `val` contains a double quote ("), the method checks if `val` also contains a space.
     * If `val` contains a space, the method returns `val` enclosed in double quotes and replaces all double quotes in `val` with three double quotes.
     * If `val` does not contain a space, the method returns `val` with all double quotes replaced with two double quotes.
     *
     * If `val` contains a space but does not contain a double quote, the method returns `val` enclosed in double quotes.
     *
     * If `val` does not contain a space or a double quote, the method returns `val` as is.
     *
     * @param {string} val - The string to be escaped.
     * @returns {string} The escaped string.
     */
    static winShellEscape(val) {
        if (val.charAt(0) === '-') {
            return val
        }
        if (~val.indexOf('"')) {
            if (~val.indexOf(' ')) {
                return `"${val.replace(DBL_QUOTE_RX, '"""')}"`
            } else {
                return val.replace(DBL_QUOTE_RX, '""')
            }
        }
        if (~val.indexOf(' ')) {
            return `"${val}"`
        }
        return val
    }

    /**
     * Parses a command string into an array of command and arguments.
     *
     * This method is used to parse a command string into an array of command and arguments.
     *
     * The method checks if the command string contains any quotes. If it doesn't, the method
     * splits the command string by spaces and returns the resulting array.
     *
     * If the command string contains quotes, the method splits the command string into an array
     * of characters and processes each character.
     *
     * If a character is a quote, the method checks if it's the same as the current quote character.
     * If it is, the method checks if the previous character is a backslash. If it is, the method
     * replaces the backslash with the quote. If it's not, the method ends the current token and
     * starts a new one.
     *
     * If a character is a space, the method checks if it's inside a quote. If it is, the method
     * adds it to the current token. If it's not, the method ends the current token and starts a new one.
     *
     * If a character is not a quote or a space, the method adds it to the current token.
     *
     * The method returns an array of tokens (command and arguments).
     *
     * @param {string} cmd - The command string to be parsed.
     * @returns {Array} An array of command and arguments.
     */
    static parseCommand(cmd) {
        if (!QUOTE_RX.test(cmd)) return cmd.split(' ')
        const chars = [...cmd]
        const lastIdx = chars.length - 1
        return chars.reduce(
            (accum, c, idx) => {
                const {tokens, token, quotes} = accum
                if (c === "'" || c === '"') {
                    if (quotes.get()) {
                        if (quotes.get() === c) {
                            if (token[token.length - 1] === '\\') {
                                token.pop()
                                token.push(c)
                            } else {
                                if (token.length) tokens.push(token.join(''))
                                token.length = quotes.clear() || 0
                            }
                        } else {
                            token.push(c)
                        }
                    } else {
                        quotes.set(undefined, c)
                    }
                } else if (c === ' ') {
                    if (quotes.get()) {
                        token.push(c)
                    } else if (token.length) {
                        tokens.push(token.join(''))
                        token.length = 0
                    }
                } else {
                    token.push(c)
                }
                if (idx === lastIdx && token.length) tokens.push(token.join(''))
                return accum
            },
            {tokens: [], token: [], quotes: new Map()}
        ).tokens
    }


    /**
     * Creates a summary of the contents of an object.
     *
     * This static method of the `CppReferenceExtension` class is used to create a summary of the contents of an object.
     * The summary is a copy of the object where all the properties whose type are Array or object are replaced with "[...]" or "{...}"
     * when the number of elements is greater than 3. Otherwise, the property is recursively replaced with its own summary.
     *
     * It takes two arguments: `obj` and `level`. `obj` is the object to be summarized. `level` is the depth level of the object properties.
     * The default value of `level` is 0.
     *
     * The method returns an object that is a summary of the input object.
     *
     * @param {Object} obj - The object to be summarized.
     * @param {number} [level=0] - The depth level of the object properties.
     * @returns {Object} An object that is a summary of the input object.
     */
    static objectSummary(obj, level = 0) {
        let summary = {}
        const maxPropertiesPerLevel = {
            0: 10,
            1: 5,
            2: 3,
            3: 2,
        }
        const maxProperties = maxPropertiesPerLevel[level] || 1

        if (Array.isArray(obj)) {
            if (obj.length > maxProperties) {
                return '[...]'
            } else {
                let arr = []
                for (const value of obj) {
                    if (typeof value === 'object') {
                        arr.push(CppReferenceExtension.objectSummary(value, level + 1))
                    } else {
                        arr.push(value)
                    }
                }
                return arr
            }
        }

        for (const key in obj) {
            let value = obj[key]
            if (Array.isArray(value)) {
                if (value.length > maxProperties) {
                    value = '[...]'
                } else {
                    value = CppReferenceExtension.objectSummary(value, level + 1)
                }
            } else if (typeof value === 'object') {
                if (Object.keys(value).length > maxProperties) {
                    value = '{...}'
                } else {
                    value = CppReferenceExtension.objectSummary(value, level + 1)
                }
            }
            summary[key] = value
        }
        return summary
    }
}

module.exports = CppReferenceExtension

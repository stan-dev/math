
$REPONAME=[io.path]::GetFileNameWithoutExtension($(git config --get remote.origin.url))
echo "REPONAME is $REPONAME"
$BOOST_SRC_FOLDER=$(git rev-parse --show-toplevel)
echo "BOOST_SRC_FOLDER is $BOOST_SRC_FOLDER"
$PARENTNAME=[io.path]::GetFileNameWithoutExtension($(git --git-dir $BOOST_SRC_FOLDER/../../.git config --get remote.origin.url))

if ( $PARENTNAME -eq "boost" )
{
    echo "Starting out inside boost-root"
    $BOOSTROOTLIBRARY="yes"
}
else
{
    echo "Not starting out inside boost-root"
    $BOOSTROOTLIBRARY="no"
}

$REPO_BRANCH=git rev-parse --abbrev-ref HEAD
echo "REPO_BRANCH is $REPO_BRANCH"

if ( $REPO_BRANCH -eq "master" )
{
    $BOOST_BRANCH="master"
}
else
{
    $BOOST_BRANCH="develop"
}

echo "BOOST_BRANCH is $BOOST_BRANCH"

echo '==================================> INSTALL'

echo "Install chocolatey"
iex ((new-object net.webclient).DownloadString('https://chocolatey.org/install.ps1'))

choco install -y rsync sed doxygen.install xsltproc docbook-bundle

if ( -Not (Get-Command java -errorAction SilentlyContinue) )
{
    choco install -y openjdk --version=17.0.1
}

if ( -Not (Get-Command python -errorAction SilentlyContinue) )
{
    choco install -y python3
}

if ( -Not (Get-Command git -errorAction SilentlyContinue) )
{
    choco install -y git
}

# Make `refreshenv` available right away, by defining the $env:ChocolateyInstall
# variable and importing the Chocolatey profile module.
# Note: Using `. $PROFILE` instead *may* work, but isn't guaranteed to.
$env:ChocolateyInstall = Convert-Path "$((Get-Command choco).Path)\..\.."
Import-Module "$env:ChocolateyInstall\helpers\chocolateyProfile.psm1"

refreshenv

# A bug fix, which may need to be developed further:
# b2 reports that the "cp" command can't be found on Windows.
# Let's add git's version of "cp" to the PATH.
$newpathitem="C:\Program Files\Git\usr\bin"
if( (Test-Path -Path $newpathitem) -and -Not ( $env:Path -like "*$newpathitem*"))
{
       $env:Path += ";$newpathitem"
}

Copy-Item "C:\Program Files\doxygen\bin\doxygen.exe" "C:\Windows\System32\doxygen.exe"

cd $BOOST_SRC_FOLDER
cd ..
if ( -Not (Test-Path -Path "tmp") )
{
    mkdir tmp
}

cd tmp

# Install saxon
if ( -Not (Test-Path -Path "C:\usr\share\java\Saxon-HE.jar") )
{
    $source = 'https://sourceforge.net/projects/saxon/files/Saxon-HE/9.9/SaxonHE9-9-1-4J.zip/download'
    $destination = 'saxonhe.zip'
    if ( Test-Path -Path $destination)
    {
        rm $destination
    }
    if ( Test-Path -Path "saxonhe")
    {
        rm Remove-Item saxonhe -Recurse -Force
    }
    Start-BitsTransfer -Source $source -Destination $destination
    Expand-Archive .\saxonhe.zip
    cd saxonhe
    if ( -Not (Test-Path -Path "C:\usr\share\java") )
    {
        mkdir "C:\usr\share\java"
    }
    cp saxon9he.jar Saxon-HE.jar
    cp Saxon-HE.jar "C:\usr\share\java\"
}

cd $BOOST_SRC_FOLDER

if ( $BOOSTROOTLIBRARY -eq "yes" )
{
    echo "updating boost-root"
    cd ../..
    git checkout $BOOST_BRANCH
    git pull
    $Env:BOOST_ROOT=Get-Location | Foreach-Object { $_.Path }
    echo "Env:BOOST_ROOT is $Env:BOOST_ROOT"
}
else
{
    cd ..
    if ( -Not (Test-Path -Path "boost-root") )
    {
        echo "cloning boost-root"
        git clone -b $BOOST_BRANCH https://github.com/boostorg/boost.git boost-root --depth 1
        cd boost-root
        $Env:BOOST_ROOT=Get-Location | Foreach-Object { $_.Path }
        echo "Env:BOOST_ROOT is $Env:BOOST_ROOT"
        if ( -Not (Test-Path -Path "libs\$REPONAME") )
        {
            mkdir libs\$REPONAME
        }
        Copy-Item -Path $BOOST_SRC_FOLDER\* -Destination libs\$REPONAME\ -Recurse -Force
    }
    else
    {
        echo "updating boost-root"
        cd boost-root
        git checkout $BOOST_BRANCH
        git pull
        $Env:BOOST_ROOT=Get-Location | Foreach-Object { $_.Path }
        echo "Env:BOOST_ROOT is $Env:BOOST_ROOT"
        if ( -Not (Test-Path -Path "libs\$REPONAME") )
        {
            mkdir libs\$REPONAME
        }
        Copy-Item -Path $BOOST_SRC_FOLDER\* -Destination libs\$REPONAME\ -Recurse -Force
    }
}

git submodule update --init libs/context
git submodule update --init tools/boostbook
git submodule update --init tools/boostdep
git submodule update --init tools/docca
git submodule update --init tools/quickbook
git submodule update --init tools/build

$matcher='\.saxonhe_jar = \$(jar\[1\]) ;$'
$replacer='.saxonhe_jar = $(jar[1]) ;  .saxonhe_jar = \"/usr/share/java/Saxon-HE.jar\" ;'
sed -i "s~$matcher~$replacer~" tools/build/src/tools/saxonhe.jam

python tools/boostdep/depinst/depinst.py ../tools/quickbook

echo "Running bootstrap.bat"
./bootstrap.bat

echo "Running ./b2 headers"
./b2 headers

$content='using doxygen : "/Program Files/doxygen/bin/doxygen.exe" ; using boostbook ; using saxonhe ;'
$filename="$Env:BOOST_ROOT\tools\build\src\user-config.jam"
[IO.File]::WriteAllLines($filename, $content)

./b2 libs/$REPONAME/doc/


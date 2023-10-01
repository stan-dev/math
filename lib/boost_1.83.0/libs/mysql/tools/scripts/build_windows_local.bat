@REM
@REM Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
@REM
@REM Distributed under the Boost Software License, Version 1.0. (See accompanying
@REM file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
@REM

SET IMAGE=build-msvc14_3
SET IMAGE_TAG=4e726c8dd2c8b4f9589a5c3ea892301e7e4a285f

SET CONTAINER="builder-%IMAGE%"
docker start %CONTAINER% || docker run -dit --name %CONTAINER% -v "%USERPROFILE%\mysql:C:\boost-mysql" "ghcr.io/anarthal-containers/%IMAGE%:%IMAGE_TAG%" || exit /b 1
docker exec %CONTAINER% python.exe "C:\boost-mysql\tools\ci.py" --source-dir=C:\boost-mysql --toolset=msvc ^
    --build-kind=cmake ^
    "--generator=Visual Studio 17 2022" ^
    --build-shared-libs=1 ^
    --clean=1 ^
    --cmake-standalone-tests=0 ^
    --cxxstd=20 ^
    --variant=debug ^
    --address-model=32 || exit /b 1

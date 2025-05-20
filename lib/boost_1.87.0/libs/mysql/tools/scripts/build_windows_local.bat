@REM
@REM Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
@REM
@REM Distributed under the Boost Software License, Version 1.0. (See accompanying
@REM file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
@REM

SET IMAGE=build-msvc14_3
SET IMAGE_TAG=latest
SET SCRIPT_PATH=%~dp0

SET CONTAINER="builder-%IMAGE%"
docker start %CONTAINER% || docker run -dit --name %CONTAINER% -v "%SCRIPT_PATH%..\..:C:\boost-mysql" "ghcr.io/anarthal-containers/%IMAGE%:%IMAGE_TAG%" || exit /b 1
docker exec %CONTAINER% python.exe "C:\boost-mysql\tools\ci\main.py" --source-dir=C:\boost-mysql b2 --toolset=msvc ^
    --cxxstd=20 ^
    --variant=debug ^
    --address-model=64 || exit /b 1

@echo off
setlocal
if not exist build-mingw mkdir build-mingw
cmake -S . -B build-mingw -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake --build build-mingw -j
echo.
echo Built: build-mingw\uagd_nsga.exe
endlocal

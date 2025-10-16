@echo off
setlocal
if not exist build mkdir build
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded
cmake --build build --config Release -j
echo.
echo Built: build\Release\uagd_nsga.exe
endlocal

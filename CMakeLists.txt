cmake_minimum_required(VERSION 3.12)
project(uagd_nsga CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
  # 静态运行库，生成更便携的 EXE（CMake>=3.15）
  if (DEFINED CMAKE_MSVC_RUNTIME_LIBRARY)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")
  else()
    add_compile_options(/MT)
  endif()
  add_compile_options(/O2 /DNOMINMAX)
else()
  add_compile_options(-O3 -march=native -funroll-loops)
endif()

add_executable(uagd_nsga uagd_nsga.cpp)

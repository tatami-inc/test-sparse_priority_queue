cmake_minimum_required(VERSION 3.14)

project(testing
    VERSION 1.0.0
    DESCRIPTION "A C++ interface for assorted matrix representations"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_executable(testing test.cpp)

include(FetchContent)

FetchContent_Declare(
  cli11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11
  GIT_TAG v2.4.1
)

FetchContent_Declare(
    nanobench
    GIT_REPOSITORY https://github.com/martinus/nanobench
    GIT_TAG v4.3.11
)

FetchContent_MakeAvailable(cli11)
FetchContent_MakeAvailable(nanobench)

target_link_libraries(testing CLI11::CLI11 nanobench::nanobench)

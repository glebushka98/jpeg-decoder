include_directories(.)
include_directories(../commons)
include_directories(../contrib)
include_directories(../contrib/gtest)
include_directories(../contrib/gmock)
include_directories(../contrib/benchmark/include)

# add_gtest: common boilerplate for creation gtest executable
function(add_gtest BINARY_NAME)
    list(REMOVE_AT ARGV 0)

    add_library(gmock ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/gmock-gtest-all.cc)

    add_executable(${BINARY_NAME}
            ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/gmock_main.cc
            ${ARGV})
    target_link_libraries(${BINARY_NAME} gmock pthread dl)
endfunction(add_gtest)

add_library(benchmark
        ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/benchmark/src/benchmark.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/benchmark/src/benchmark_register.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/benchmark/src/colorprint.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/benchmark/src/commandlineflags.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/benchmark/src/complexity.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/benchmark/src/console_reporter.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/benchmark/src/counter.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/benchmark/src/csv_reporter.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/benchmark/src/json_reporter.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/benchmark/src/reporter.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/benchmark/src/sleep.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/benchmark/src/statistics.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/benchmark/src/string_util.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/benchmark/src/sysinfo.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/benchmark/src/timers.cc)
target_include_directories(benchmark
        PRIVATE ../contrib/benchmark/src)
target_compile_definitions(benchmark
        PRIVATE HAVE_POSIX_REGEX)

# add_benchmark: common boilerplate for benchmark executable
function(add_benchmark BINARY_NAME)
    list(REMOVE_AT ARGV 0)

    add_executable(${BINARY_NAME} ${ARGV})
    target_link_libraries(${BINARY_NAME} benchmark pthread)
endfunction(add_benchmark)

set(CXX_STANDARD 17)
if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.1)
    set(CXX_STANDARD 14)
endif()

if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
    set(CXX_STANDARD 11)
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND APPLE)
    set(CXX_STANDARD 14)
endif()

if (CMAKE_COMPILER_IS_GNUCC OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${CXX_STANDARD} -Wall -Wextra -Wpedantic -g")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -g")
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fno-omit-frame-pointer -Wall -Wextra -Wpedantic -g -O0")

set(CMAKE_CXX_FLAGS_ASAN "-g -fsanitize=address,undefined -fno-sanitize-recover=all"
        CACHE STRING "Compiler flags in asan build"
        FORCE)

set(CMAKE_CXX_FLAGS_TSAN "-g -fsanitize=thread -fno-sanitize-recover=all"
        CACHE STRING "Compiler flags in tsan build"
        FORCE)

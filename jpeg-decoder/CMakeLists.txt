cmake_minimum_required(VERSION 3.0)
project(jpeg-decoder)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIRS})

find_package(JPEG REQUIRED)
include_directories(${JPEG_INCLUDE_DIR})

find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})

include(../common.cmake)

if (TEST_SOLUTION)
  include_directories(../private/jpeg-decoder)

  add_library(decoder-lib
          ../private/jpeg-decoder/decoder.cpp
          ../private/jpeg-decoder/huffman.cpp)
else()
  # add your libraries here
  add_library(decoder-lib decoder.cpp huffman.cpp decoder_utils.cpp)
endif()


add_executable(test_baseline
        test_baseline.cpp
        ../contrib/catch_main.cpp)

target_link_libraries(test_baseline decoder-lib ${FFTW_LIBRARIES} ${PNG_LIBRARY} ${JPEG_LIBRARIES})

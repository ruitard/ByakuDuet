cmake_minimum_required(VERSION 3.0.0)
project(ByakuDuet VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-fsanitize=address)
add_link_options(-fsanitize=address)

find_package(OpenCV CONFIG REQUIRED)

add_subdirectory(opencv_dxt)
add_subdirectory(src)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

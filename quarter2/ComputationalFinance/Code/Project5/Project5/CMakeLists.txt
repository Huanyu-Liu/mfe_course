cmake_minimum_required(VERSION 3.0.0)
project(project5 VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

add_executable(project5 main.cpp option.h option.cpp laguerre.h laguerre.cpp matrix_cal.h matrix_cal.cpp random_generator.h random_generator.cpp stats.h stats.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

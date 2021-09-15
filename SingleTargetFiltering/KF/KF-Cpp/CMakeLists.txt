cmake_minimum_required(VERSION 3.5)
project(Kalman_Filter)

set(CMAKE_CXX_STANDARD 11)

find_package(Boost REQUIRED)
#include_directories(${Boost_INCLUDE_DIR})

set(SOURCE_FILES KalmanFilter.cpp KalmanFilter_demo.cpp)
add_executable(KF_demo ${SOURCE_FILES})
target_link_libraries(KF_demo ${BOOST_LIBRARIES})


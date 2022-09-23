cmake_minimum_required(VERSION 3.23)
project(VDFs)

set(CMAKE_CXX_STANDARD 20)

set(SOURCES main.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE NTL/include)

add_subdirectory(NTL/src)
target_link_libraries(${PROJECT_NAME} ntl_lib)

set(LIBRARY_NAME "ntl_lib")

file(GLOB SOURCES "*.cpp")

add_library(${LIBRARY_NAME} ${SOURCES})

target_include_directories(${LIBRARY_NAME} PRIVATE ../include)

cmake_minimum_required(VERSION 3.9)
project(lammpsreader Fortran)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fcheck=all")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -march=native")
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

add_library(lammpsreader
    STATIC
    reader.f90
)

add_subdirectory(coretran/src)

target_link_libraries(lammpsreader coretran)

target_include_directories(lammpsreader PUBLIC ${PROJECT_BINARY_DIR}/coretran/include)

install(TARGETS lammpsreader DESTINATION lib)
install(DIRECTORY ${PROJECT_BINARY_DIR}/include/ DESTINATION include)

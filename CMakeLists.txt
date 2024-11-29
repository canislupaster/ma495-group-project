cmake_minimum_required(VERSION 3.29)
project(mathmodel C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_COLOR_DIAGNOSTICS ON)
add_compile_options(-Wall -Wextra -Wpedantic -Wextra -Wshadow -Wno-sign-compare)

if (CMAKE_BUILD_TYPE MATCHES Release)
    add_compile_options(-Ofast)
endif()

if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_compile_definitions(BUILD_DEBUG)
    add_compile_options(-fsanitize=address,undefined)
    add_link_options(-fsanitize=address,undefined)
endif ()

add_library(gsw ./GSW-C/gsw_oceanographic_toolbox.c ./GSW-C/gsw_saar.c)

find_package(HDF5 REQUIRED)
find_package(GeographicLib REQUIRED)

include(FetchContent)
FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
FetchContent_MakeAvailable(json)

add_executable(main main.cpp)
target_link_libraries(main HDF5::HDF5 ${GeographicLib_LIBRARIES} gsw nlohmann_json::nlohmann_json)
target_include_directories(main PRIVATE HighFive/include ${GeographicLib_INCLUDE_DIRS} ./GSW-C)
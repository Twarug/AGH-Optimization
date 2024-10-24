cmake_minimum_required(VERSION 3.29)
project(Opt)

set(CMAKE_CXX_STANDARD 26)

file(GLOB_RECURSE SOURCE CONFIGURE_DEPENDS src/**.cpp)
file(GLOB_RECURSE HEADER CONFIGURE_DEPENDS src/**.h)

add_executable(Opt ${SOURCE})

target_sources(Opt PRIVATE ${ARCHIMEDES_SOURCE})
target_include_directories(Opt PRIVATE src)


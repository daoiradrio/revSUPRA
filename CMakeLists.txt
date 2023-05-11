cmake_minimum_required(VERSION 3.16)
project(SUPRA)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Debug)

add_compile_options(-fsanitize=address)
add_link_options(-fsanitize=address)

include_directories(include)

file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable(supra main_files/main.cpp ${SRC_FILES})
# add more executables?

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

enable_testing()

# add tests like:
# add_executable(queue_tests tests/queue_tests.cpp ${SRC_FILES})
# target_link_libraries(queue_tests GTest::gtest_main)
# gtest_discover_tests(queue_tests)

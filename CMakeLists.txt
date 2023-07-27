cmake_minimum_required(VERSION 3.26)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

project(
    "PolynomialsCXX" 
    VERSION 1.0
    DESCRIPTION ""
    HOMEPAGE_URL ""
    )

option(POLYNOMIALS_CLONE_SUBMODULE_EIGEN "SET TO TRUE TO CLONE THE EIGEN REPOSITORY INSIDE THE PROJECT" OFF)

if (POLYNOMIALS_CLONE_SUBMODULE_EIGEN)
    message(STATUS "> Cloning Eigen from remote...")
    execute_process(
        COMMAND git submodule update --init eigen
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    message(STATUS "> Finished cloning Eigen from remote...")
    add_subdirectory(eigen)
    set_property(TARGET eigen PROPERTY FOLDER "Eigen3")
else ()
    find_package (Eigen3 3.3 REQUIRED NO_MODULE)
endif ()

if(NOT TARGET Eigen3::Eigen)
    message(FATAL_ERROR "Eigen3 library not found, install it before configuring the project or consider setting POLYNOMIALS_CLONE_SUBMODULE_EIGEN to true")
endif()

add_library(polynomials STATIC polynomial.cpp polynomial.hpp)
target_include_directories(polynomials PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(polynomials PRIVATE Eigen3::Eigen)
set_property(TARGET polynomials PROPERTY FOLDER "Polynomials")

option(POLYNOMIALS_WITH_TESTS "Compile tests" ON)

if (POLYNOMIALS_WITH_TESTS)
    enable_testing ()
    
    add_executable(test_polynomial test_polynomial.cpp)
    target_link_libraries(test_polynomial PUBLIC polynomials)
    target_include_directories(test_polynomial PRIVATE $(CMAKE_CURRENT_SOURCE_DIR))
    add_test(NAME test_polynomial COMMAND test_polynomial)

    set_property(TARGET test_polynomial PROPERTY FOLDER "Polynomials")
endif ()

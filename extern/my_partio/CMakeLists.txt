cmake_minimum_required(VERSION 3.21)

file(GLOB io_cpp "io/*.cpp")
file(GLOB core_cpp "core/*.cpp")

add_library(mypartio  ${io_cpp} ${core_cpp})
set_target_properties(mypartio PROPERTIES
    OUTPUT_NAME mypartio
    POSITION_INDEPENDENT_CODE ON
    CXX_STANDARD 17
)

target_compile_features(mypartio PUBLIC cxx_std_17)

IF(MSVC)
  target_compile_definitions(mypartio 
  PRIVATE 
    PARTIO_WIN32=1
    _USE_MATH_DEFINES=1
    )
    target_compile_options(mypartio PRIVATE /W0)

ELSE(MSVC)
  target_compile_options(mypartio 
  PRIVATE  
    -Wextra 
    -Wno-unused-parameter
    -fPIC
    -w)

ENDIF(MSVC)

target_compile_definitions(mypartio  PUBLIC PARTIO_USE_ZLIB=1)


target_include_directories(mypartio PUBLIC ${CMAKE_CURRENT_LIST_DIR})

# find_package(ZLIB REQUIRED)
# if (ZLIB_FOUND)
#     target_link_libraries(mypartio PUBLIC ZLIB::ZLIB)
# endif()

target_include_directories(mypartio PUBLIC ../zlib/src )


option(TEST_MYPARTIO "test my partio" OFF)
if(TEST_MYPARTIO)
# test_mypartio
    add_executable(test_mypartio my_test/test_mypartio.cpp)
    target_compile_definitions(test_mypartio PRIVATE  PARTIO_USE_ZLIB)
    target_link_libraries(test_mypartio mypartio )
    add_test(
      NAME test_mypartio
      COMMAND $<TARGET_FILE:test_mypartio>
    )

# test_readwrite_particles
    add_executable(test_readwrite_particles my_test/test_readwrite_particles.cpp)
    target_compile_definitions(test_readwrite_particles PRIVATE  PARTIO_USE_ZLIB)
    target_link_libraries(test_readwrite_particles mypartio )

# test_myio
    add_executable(test_myio my_test/test_myio.cpp)
    target_compile_definitions(test_myio PRIVATE  PARTIO_USE_ZLIB)
    target_link_libraries(test_myio mypartio )

# test_readWriteVector.cpp
  add_executable(test_readWriteVector my_test/test_readWriteVector.cpp)
  target_compile_definitions(test_readWriteVector PRIVATE PARTIO_USE_ZLIB)
  target_include_directories(test_myio PRIVATE  ${EIGEN3_INCLUDE_DIR})
  target_link_libraries(test_readWriteVector mypartio)
endif()

################################################################################
# general CMake and project setup
################################################################################
cmake_minimum_required(VERSION 3.23)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(SPlisHSPlasH)

# MYADD 关闭警告C4819: 该文件包含不能在当前代码页(936)中表示的字符。请将该文件保存为 Unicode 格式以防止数据丢失
#  关闭警告C4267: “参数”: 从“size_t”转换到“int”，可能丢失数据
if(MSVC)
	add_compile_options(/wd4819 /wd4267)
endif()

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
include(Common)
if (NOT WIN32)
	message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
endif()

set(ExternalInstallDir "${CMAKE_BINARY_DIR}/extern/install" CACHE INTERNAL "")
set(EXT_CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE INTERNAL "")
if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
	set(EXT_CMAKE_BUILD_TYPE "Release" CACHE INTERNAL "")
endif()

set(TOPLEVEL_INCLUDE_DIR ${PROJECT_SOURCE_DIR})

option(USE_IMGUI "Use imgui instead of AntTweakBar"	ON)


# MYADD： 增加Just My Code 方便调试
if(MSVC)
	# if($<CONFIG:Debug>)
		message(STATUS "Just My Code")
		set(VS_JUST_MY_CODE_DEBUGGING ON)
	# 	add_compile_options( /JMC )
	# endif()
endif()

################################################################################
# foreign external libraries
################################################################################
set(EXTERN_DIR "${PROJECT_SOURCE_DIR}/extern" CACHE INTERNAL "")

# MYADD： 防止target发生冲突
include(CMake/target_clash_guard.cmake)
add_subdirectory(extern/zlib)
add_subdirectory(extern/partio)
add_subdirectory(extern/my_partio) #MYADD
# add_subdirectory(extern/md5)
target_clash_guard(md5 extern/md5)
add_subdirectory(extern/tinyexpr)
if (NOT SPH_LIBS_ONLY)
	# add_subdirectory(extern/AntTweakBar)
	# add_subdirectory(extern/glfw)
	target_clash_guard(AntTweakBar extern/AntTweakBar)
	target_clash_guard(glfw extern/glfw)
	add_subdirectory(extern/imgui)
	if (USE_PYTHON_BINDINGS)
		add_subdirectory(extern/pybind)
	endif ()
endif()

## Eigen3 is used by most of the libraries that follow
# include(FindEigen3)
set(Eigen3_DIR "${PROJECT_SOURCE_DIR}/extern/eigen")
find_package(Eigen3 REQUIRED)
add_definitions(-DEIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT)

################################################################################
# own external libraries
################################################################################
include(SetUpExternalProjects)

################################################################################
# internal libraries
################################################################################
add_subdirectory(SPlisHSPlasH)
add_subdirectory(Utilities)

################################################################################
# executables
################################################################################
if (NOT SPH_LIBS_ONLY)
	include(DataCopyTargets)
	add_subdirectory(Tools)
	add_subdirectory(Simulator)

	# MYADD
	option(ENABLE_TESTS "enable the tests" OFF)
	if(ENABLE_TESTS)
		enable_testing()
		add_subdirectory(Tests)
	endif()

	if (USE_PYTHON_BINDINGS)
		add_subdirectory(pySPlisHSPlasH)
	endif ()
endif()

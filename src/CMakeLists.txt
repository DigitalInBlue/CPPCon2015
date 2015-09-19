#
# Top level makefile for "Demystifying Floating Point", CPPCon 2015
#
# Copyright 2015 John Farrier 
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#

#
# Cmake Configuration
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project(CPPCon2015)

include(CheckFunctionExists)
include(CheckCXXSourceCompiles)
include(CheckIncludeFile)

#
# User Options
#

if(DFP_COMPILE_DYNAMIC_LIBRARIES)
	SET(DFP_USER_DEFINED_SHARED_OR_STATIC "SHARED")
else()
	SET(DFP_USER_DEFINED_SHARED_OR_STATIC "STATIC")
endif()

set(GTEST_DIRECTORY gtest-1.7.0)
set(GTEST_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${GTEST_DIRECTORY}/include)
set(GTEST_LIBRARY gtest)
set(GTEST_MAIN_LIBRARY gtest_main)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${GTEST_DIRECTORY}/include ${CMAKE_CURRENT_SOURCE_DIR}/${GTEST_DIRECTORY}/src)
ADD_SUBDIRECTORY(${GTEST_DIRECTORY})

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#
# Compiler Settings
#

if(MSVC) 
  # Force to always compile with warning level 4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]") 
	string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}") 
  else() 
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4") 
  endif() 

  set(CPPCON2015_FLOATING_POINT_MODE "/fp:precise" CACHE STRING "Compiler floating point flag.")
  set_property(CACHE CPPCON2015_FLOATING_POINT_MODE PROPERTY STRINGS "/fp:precise" "/fp:fast" "/fp:strict")

  add_definitions("/D_CRT_SECURE_NO_WARNINGS /wd4251 /nologo")
  add_definitions(${CPPCON2015_FLOATING_POINT_MODE})
endif() 
 
if(BORLAND) 
  set(CMAKE_CXX_FLAGS "" CACHE STRING "Compiler Flags for All Builds" FORCE) 
  set(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING "Compiler Flags for Debug Builds" FORCE)  
  set(CMAKE_CXX_FLAGS_RELEASE  "" CACHE STRING "Compiler Flags for Release Builds" FORCE) 
endif() 

message (STATUS "SYSTEM: ${CMAKE_SYSTEM_NAME}")
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++") 
	set(CMAKE_COMPILER_IS_GNUCXX 1)
endif()
 
if(UNIX)
	if(CMAKE_COMPILER_IS_GNUCXX) 
	  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -std=gnu++11 -Wall -O3" CACHE STRING "Compiler Flags for All Builds" FORCE) 
	  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pedantic -O0 -g" CACHE STRING "Compiler Flags for Debug Builds" FORCE) 
	  set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "Compiler Flags for Release Builds" FORCE) 
	else()
	  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -std=gnu++11 -Wall -fPIC -O3" CACHE STRING "Compiler Flags for All Builds" FORCE) 
	  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pedantic -O0 -g" CACHE STRING "Compiler Flags for Debug Builds" FORCE) 
	  set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "Compiler Flags for Release Builds" FORCE) 
	endif() 
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#
# Build and Install Settings
#

set(CMAKE_DEBUG_POSTFIX          "d" CACHE STRING "add a postfix, usually d on windows")
set(CMAKE_RELEASE_POSTFIX        ""  CACHE STRING "add a postfix, usually empty on windows")
set(CMAKE_RELWITHDEBINFO_POSTFIX ""  CACHE STRING "add a postfix, usually empty on windows")
set(CMAKE_MINSIZEREL_POSTFIX     ""  CACHE STRING "add a postfix, usually empty on windows")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

message("CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")

# --------------------------------------------------------------------------- 
# Google Test Application
# --------------------------------------------------------------------------- 

set(PROJECT_NAME CPPCon2015)

set(GTEST_DIRECTORY gtest-1.7.0)
set(GTEST_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${GTEST_DIRECTORY}/include)
set(GTEST_LIBRARY gtest)
set(GTEST_MAIN_LIBRARY gtest_main)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${GTEST_DIRECTORY}/include ${CMAKE_CURRENT_SOURCE_DIR}/${GTEST_DIRECTORY}/src)

add_executable(${PROJECT_NAME} 
	CPPCon2015.test.cpp
)

# VS2012 doesn't support true variadic templates
if(MSVC) 
	add_definitions( /D _VARIADIC_MAX=10 )
endif()

include_directories(${HEADER_PATH})
include_directories(${GTEST_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})


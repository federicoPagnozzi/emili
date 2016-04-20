IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Debug)
ENDIF()

project(EMILI)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)

aux_source_directory(./examtt EXAMTT_SRC_LIST)
aux_source_directory(./SA SA_SRC_LIST)
aux_source_directory(./QAP QAP_SRC_LIST)
aux_source_directory(./pfsp PFSP_SRC_LIST)

option(USE_SSE "Use SSE compile path if available" OFF)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-sign-compare -Wno-unused -Wno-parentheses")

#SSE
if(USE_SSE)
    message(STATUS "Using sse compile path")
    add_definitions(-DENABLE_SSE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse")
endif(USE_SSE)

#C++ 11 detection
IF("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)

    IF(NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
       add_definitions(-DNOC11)
       message(WARNING "Your gcc does not support c++11, switched to c++0x support.")
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
    ELSE()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    ENDIF()
ELSEIF("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ELSE()
    add_definitions(-DNOC11)
    message(WARNING, "Your gcc does not support c++11, switched to c++0x support.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    # message(FATAL_ERROR "Your C++ compiler does not support C++11.")
ENDIF()

SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")

add_executable(${PROJECT_NAME} emilibase.cpp ${EXAMTT_SRC_LIST} ${PFSP_SRC_LIST} generalParser.cpp main.cpp ${SA_SRC_LIST} ${QAP_SRC_LIST})


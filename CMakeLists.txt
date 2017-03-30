project(emili)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
aux_source_directory(./pfsp SRC_LIST)
aux_source_directory(./vrp SRC_LIST)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -std=c++11")
option(USE_SSE "Use SSE compile path if available" OFF)
option(USE_NEW_MAIN "Use new loading paradigm" OFF)
option(WITH_STATS "Prints run informations" OFF)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -O3")


#ADD COMMIT HASH TO SOURCE CODE
add_custom_target(GNUM ALL git rev-parse HEAD | xargs -I cnum sed -ie s/.*"define GIT_COMMIT_NUMBER".*/"\#define GIT_COMMIT_NUMBER \\"cnum\\""/g generalParser.cpp WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_custom_target(CGNUM ALL rm -f generalParser.cppe WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

#SSE
if(USE_SSE)
	message(STATUS "Using sse compile path") 
	add_definitions(-DENABLE_SSE)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse")
endif(USE_SSE)

#Solution stats
if(WITH_STATS)
    message(STATUS "Solution stats enabled")
    add_definitions(-DWITH_STATS)
endif(WITH_STATS)

if(USE_NEW_MAIN)
    message(STATUS "Component based algorithm loading enabled")
    add_definitions(-DMAIN_NEW)
endif(USE_NEW_MAIN)

#C++ 11 detection
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process( COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
           add_definitions(-DNOC11)
           message(WARNING "Your gcc does not support c++11, switched to c++0x support.")
           set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
    else ()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else ()
    add_definitions(-DNOC11)
    message(WARNING, "Your gcc does not support c++11, switched to c++0x support.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    #message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

add_executable(${PROJECT_NAME} ${SRC_LIST})
add_dependencies(${PROJECT_NAME} GNUM)
add_dependencies(CGNUM GNUM)

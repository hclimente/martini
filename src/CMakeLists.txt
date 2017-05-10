cmake_minimum_required(VERSION 3.6)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 11)

file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/CEasyGWAS/*.cc ${PROJECT_SOURCE_DIR}/src/*.h
		${PROJECT_SOURCE_DIR}/src/CEasyGWAS/*/*.cc ${PROJECT_SOURCE_DIR}/src/*/*.h
		${PROJECT_SOURCE_DIR}/src/External/eigen/*.cc ${PROJECT_SOURCE_DIR}/src/External/eigen/*.h)

if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
	SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g" )
endif (CMAKE_BUILD_TYPE EQUAL "DEBUG")

include_directories(${PROJECT_SOURCE_DIR}/src/
		${PROJECT_SOURCE_DIR}/External ${PROJECT_SOURCE_DIR}/External/eigen
		${PROJECT_SOURCE_DIR}/External/maxflow ${PROJECT_SOURCE_DIR}/External/Cephes)

add_library(CEasyGWAS STATIC ${SOURCES})
target_link_libraries (CEasyGWAS cephes maxflow)

add_subdirectory (${PROJECT_SOURCE_DIR}/External/Cephes)
add_subdirectory (${PROJECT_SOURCE_DIR}/External/maxflow)

install(TARGETS CEasyGWAS
		RUNTIME DESTINATION bin
		ARCHIVE DESTINATION lib)

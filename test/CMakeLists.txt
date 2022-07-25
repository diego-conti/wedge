cmake_minimum_required(VERSION 3.10)
project(TestWedge)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_definitions(IN_GINAC) 
add_compile_options(-g -O0)
link_libraries(ginac wedge cln)
link_directories(${GINAC_DIR}/ginac/.libs ${CMAKE_BINARY_DIR}/wedge)
include_directories(${GINAC_DIR} ${GINAC_DIR}/ginac ${CMAKE_SOURCE_DIR}/wedge ${CMAKE_BINARY_DIR}/wedge ${CMAKE_SOURCE_DIR}/cxxtest)

#set (CXXTESTFLAGS --error-printer --no-eh)
set (CXXTESTFLAGS --error-printer)

enable_testing()
foreach(test  algebra base connections convenience coordinates functions lambda liegroups manifolds poly representations structures wedgealgebraic)
	set (runner run${test}.cpp)
	set (testdefinition ${test}.h)
	add_executable(${test} ${runner} src/${testdefinition} src/test.h)
	add_custom_command(
	  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${runner}
	  COMMAND ${CXXTEST_DIR}/bin/cxxtestgen ${CXXTESTFLAGS} -o ${runner} ${PROJECT_SOURCE_DIR}/src/${testdefinition}
	  DEPENDS src/${testdefinition} src/test.h 
	 )
	add_test(NAME ${test} COMMAND ${test}) 
endforeach()




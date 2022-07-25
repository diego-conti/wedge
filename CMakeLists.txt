include(FetchContent)
cmake_minimum_required(VERSION 3.10)
project(Wedge VERSION 4.0)
configure_file(wedge/wedge/wedgeconfig.h.in wedge/wedgeconfig.h)
set(GINAC_DIR ${CMAKE_SOURCE_DIR}/ginac-1.8.3)
set(CXXTEST_DIR ${CMAKE_SOURCE_DIR}/cxxtest) 

FetchContent_Declare(ginac
	URL https://www.ginac.de/ginac-1.8.3.tar.bz2 
	SOURCE_DIR ${GINAC_DIR}
	PATCH_COMMAND patch -i ${CMAKE_SOURCE_DIR}/function.hppy.patch ${GINAC_DIR}/ginac/function.hppy     
	COMMAND patch -i ${CMAKE_SOURCE_DIR}/function.cppy.patch ${GINAC_DIR}/ginac/function.cppy   
	COMMAND patch -i ${CMAKE_SOURCE_DIR}/operators.cpp.patch ${GINAC_DIR}/ginac/operators.cpp   
	COMMAND patch -i ${CMAKE_SOURCE_DIR}/print.h.patch ${GINAC_DIR}/ginac/print.h       
	COMMAND patch -i ${CMAKE_SOURCE_DIR}/ncmul.cpp.patch ${GINAC_DIR}/ginac/ncmul.cpp  
	COMMAND patch -i ${CMAKE_SOURCE_DIR}/CMakeLists.txt.patch ${GINAC_DIR}/CMakeLists.txt  
)
FetchContent_Populate(ginac)

FetchContent_Declare(cxxtest
	GIT_REPOSITORY  https://github.com/CxxTest/cxxtest
	GIT_TAG 191adddb3876ab389c0c856e1c03874bf70f8ee4
	GIT_SHALLOW true 
	SOURCE_DIR ${CXXTEST_DIR}
	CONFIGURE_COMMAND ""
	TEST_COMMAND ""
	BUILD_COMMAND ""
	INSTALL_COMMAND ""
)
FetchContent_Populate(cxxtest)



add_subdirectory(${GINAC_DIR} ${GINAC_DIR}) #in-source build 
add_subdirectory(wedge)
add_subdirectory(test)
add_subdirectory(examples)
#install(TARGETS ginac PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#install(DIRECTORY ${GINAC_DIR}/ginac/ DESTINATION include/wedge/ginac FILES_MATCHING PATTERN "*.h*")

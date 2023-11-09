#Wedge

A [GiNaC](http://www.ginac.de)-based C++ library for symbolic computations in differential geometry.

###Requirements

You need to have installed [cmake](https://cmake.org/), [gmp](https://gmplib.org/), [cln](https://www.ginac.de/CLN) and [CoCoa](https://cocoa.dima.unige.it/cocoa/cocoalib/).
*Wedge* also uses [cxxtest](http://cxxtest.com) and [GiNaC](http://www.ginac.de), which are downloaded by `cmake` at configure time.

*Wedge* 0.4.2 has been developed and tested on Ubuntu Linux 22.10 with:

* `Cocoalib` 0.99800
* `Ginac` 1.8.4
* `gcc` 12.2
* `CMake` 3.24.2

###How to build/install

	mkdir build
	cd build
	cmake ..
	cmake --build .
	cmake --install . --prefix /the/directory/where/you/want/wedge
Specifying a prefix is optional. If none is indicated, *Wedge* and a patched version of GiNaC will be installed in the default path for your system, which probably requires administrator rights.

###How to test
After building, run the following:

	cd build
	ctest --test-dir test
This command will run some automated tests to verify that everything is working.
###Examples

Some examples are built together with *Wedge*. The source code can be found in the directory `examples` and the executables are built in `build/examples`.

###Using wedge

In order to compile a program using *Wedge*, create a file called `CMakeLists.txt` in the folder where the code is contained, and paste in the following content.

	cmake_minimum_required(VERSION 3.10)
	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_STANDARD_REQUIRED True)
	add_compile_options(-g -Wctor-dtor-privacy -Wreorder -Wold-style-cast -Wsign-promo -Wchar-subscripts -Winit-self -Wmissing-braces -Wparentheses -Wreturn-type -Wswitch -Wtrigraphs -	Wextra -Wno-sign-compare -Wno-narrowing -Wno-attributes)
	project(MyWedgeProject)
	set(SRC main.cpp)
	add_executable(mywedgeproject ${SRC})
	target_link_libraries(mywedgeproject PUBLIC wedge ginac cocoa gmp)
	target_link_directories(mywedgeproject PUBLIC $ENV{WEDGE_PATH}/lib)
	target_include_directories(mywedgeproject PUBLIC $ENV{WEDGE_PATH}/include)

Replace `main.cpp` with the name(s) of the C++ source files that need to be compiled. Change the compile options as per your preferences.

In your `main.cpp`, include *Wedge* with

	#include "wedge/wedge.h"
	
Now your program will compile by running the following:

	mkdir build
	cd build
	cmake ..
	cmake --build .

If you did not install *Wedge* in the default path, before running `cmake` you should set the environment variable `WEDGE_PATH` to point to that path, e.g.

	export WEDGE_PATH=/the/directory/where/you/want/wedge
	

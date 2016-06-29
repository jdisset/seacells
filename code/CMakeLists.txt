cmake_minimum_required(VERSION 2.8)
project(Alife15)

set(CMAKE_CXX_FLAGS "-O3 -g -std=c++14 -Wall -Wextra -Wundef -Wshadow -pedantic")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

file(GLOB SRC
	"src/*.h"
	"src/*.hpp"
	"src/core/*.h"
	"src/viewer/*.hpp"
	"src/core/*.hpp"
	"src/core/*.cpp"
	"src/GRN/*.h"
	"src/GRN/*.hpp"
	"src/GRN/*.cpp"
	"src/external/*.hpp"
	)

add_executable(console ${SRC} src/mainconsole.cpp)
target_link_libraries(console mecacell)

#add_definitions(-DCLUSTER)
add_definitions(-DOMP)
find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
add_executable(evo ${SRC} src/mainevo.cpp)
target_link_libraries(evo mecacell)
target_link_libraries(evo ${MPI_CXX_LIBRARIES})



set(CMAKE_AUTOMOC ON)
find_package(Qt5Core)
find_package(Qt5Widgets)
find_package(Qt5Gui)
find_package(Qt5Quick)
find_package(Qt5OpenGL)
add_executable(viewer ${SRC} src/mainviewer.cpp)
target_link_libraries(viewer mecacell)
target_link_libraries(viewer mecacellviewer)
qt5_use_modules(viewer Quick Core Gui OpenGL)


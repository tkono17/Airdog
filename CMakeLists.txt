#-----------------------------------------------------------------------
# CMakeLists.txt
#-----------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8)
set(PROJECT_NAME Airdog)

#project(${PROJECT_NAME})
#set(${PROJECT_NAME} VERSION MAJOR 1)
#set(${PROJECT_NAME} VERSION MINOR 0)

set(DESTINATION $ENV{SWDIR})
set(CMAKE_INSTALL_PREFIX ${DESTINATION})

add_subdirectory(AirLattice)


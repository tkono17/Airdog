cmake_minimum_required(VERSION 2.8)
set(PROJECT_NAME AirLattice)

project(${PROJECT_NAME})
set(${PROJECT_NAME} VERSION MAJOR 1)
set(${PROJECT_NAME} VERSION MINOR 0)
set(DESTINATION $ENV{SWDIR})

set(CMAKE_INSTALL_PREFIX ${DESTINATION})

#message(SWDIR = $ENV{SWDIR})
#message(DESTINATION = ${DESTINATION})
#message(CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX})

set(CMAKE_MODULE_PATH /nfs/space1/tkohno/sw/etc/cmake ${CMAKE_MODULE_PATH})

configure_file(
  "${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/${PROJECT_NAME}Config.h.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}/${PROJECT_NAME}Config.h"
)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

#####-------------------------------------------------------------
##### Package dependencies
#####-------------------------------------------------------------
# ROOT stuffs
set(ROOT_COMPONENTS_TO_USE Core RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread)
find_package(ROOT 
  COMPONENTS ${ROOT_COMPONENTS_TO_USE} REQUIRED)
include(${ROOT_USE_FILE})
list(APPEND LIBS ${ROOT_LIBRARIES})

file(GLOB LIB_SRCS src/[a-zA-Z]*.cxx)
file(GLOB EXE_SRCS src/main/*.cxx)
# file(GLOB DICT_HEADERS ${PROJECT_NAME}/[a-zA-Z].hxx)
set(DICT_HEADERS )

#####-------------------------------------------------------------
##### ROOT dictionary
#####-------------------------------------------------------------
if (DICT_HEADERS)
  include_directories(${PROJECT_SOURCE_DIR}/${PROJECT_NAME})
  ROOT_GENERATE_DICTIONARY(${PROJECT_NAME}Dict 
    ${DICT_HEADERS}
    LINKDEF ${PROJECT_NAME}/${PROJECT_NAME}LinkDef.hxx
  )
  set(root_dict_src ${PROJECT_NAME}Dict)
else()
  set(root_dict_src )
endif()

#####-------------------------------------------------------------
##### Shared library
#####-------------------------------------------------------------
add_library(${PROJECT_NAME} SHARED 
  ${LIB_SRCS} ${root_dict_src}
  )
target_link_libraries(${PROJECT_NAME} ${LIBS})

#####-------------------------------------------------------------
##### Executables
#####-------------------------------------------------------------
set(EXES )
foreach(progsrc ${EXE_SRCS})
  string(REGEX REPLACE ".*src/main/(.*).cxx" "\\1.exe" progexe ${progsrc})
  list(APPEND EXES ${progexe})
  add_executable(${progexe} ${progsrc})
  target_link_libraries(${progexe} ${PROJECT_NAME} ${LIBS})
endforeach()

#####-------------------------------------------------------------
##### Installation
#####-------------------------------------------------------------
install(DIRECTORY ${PROJECT_NAME}
  DESTINATION include
  FILES_MATCHING PATTERN "*.hxx"
  )
install(TARGETS ${PROJECT_NAME} DESTINATION lib)
foreach(progexe ${EXES})
  install(TARGETS ${progexe} DESTINATION bin)
endforeach()

#install(FILES share/${PROJECT_NAME}Config.cmake DESTINATION share/cmake)


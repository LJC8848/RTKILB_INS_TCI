cmake_minimum_required(VERSION 3.20)   # 版本按你的需要写
project(MyProjectName)
include_directories(.)
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
set(RUNTIME_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

if(WIN32)
  add_definitions(-DWIN_DLL)
else()
  add_definitions(-DDLL)
endif()

if (GCC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O3 -ansi -pedantic")
    set(CMAKE_C_FLAGS "-Wno-unused-but-set-variable -Wno-format-overflow -Wno-unused-result -Wpointer-to-int-cast")
endif()

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    if (GCC)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -fsee -fomit-frame-pointer -fno-signed-zeros -fno-math-errno -funroll-loops")
    endif()
endif()

if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

aux_source_directory(. DIR_SRCS_RTKLIB)
set(DIR_SRCS_RTKLIB_RCV
  ${DIR_SRCS_RTKLIB_RCV}
  rcv/binex.c
  rcv/crescent.c
  rcv/javad.c
  rcv/novatel.c
  rcv/nvs.c
  rcv/rt17.c
  rcv/septentrio.c
  rcv/skytraq.c
  rcv/swiftnav.c
  rcv/ublox.c
  rcv/unicore.c
)
list(APPEND DIR_SRCS ${DIR_SRCS_RTKLIB} ${DIR_SRCS_RTKLIB_RCV})

# loosely-coupled GNSS/SINS, port of GREAT-MSF's LCI: objects of this library
# are merged into rtklib below (postpos.c calls ins_feed() from procpos()).
add_subdirectory(ins)

add_library(rtklib SHARED ${DIR_SRCS} rtklib.h)
target_sources(rtklib PRIVATE $<TARGET_OBJECTS:t_ins>)
set_property(TARGET rtklib PROPERTY C_STANDARD 99)
set_target_properties(rtklib PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)


target_include_directories(rtklib
    PUBLIC ${PROJECT_SOURCE_DIR}/include
)

if(WIN32)
  target_link_libraries(rtklib wsock32 ws2_32 winmm)
  add_definitions(-DWIN_DLL)
else()
  target_link_libraries(rtklib m pthread)
  add_definitions(-DDLL)
endif()
install(TARGETS rtklib RUNTIME DESTINATION lib)

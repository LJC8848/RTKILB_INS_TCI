# src/ins/CMakeLists.txt
#
# Loosely-coupled GNSS/SINS (plain C, port of GREAT-MSF's LCI).  Built as an
# OBJECT library and merged into the rtklib DLL by src/CMakeLists.txt, because
# the post-processing hook lives in postpos.c, whose helpers (procpos, rtkpos,
# the navs/obs statics) are static to that library - a separate DLL could not
# reach them.
#
# WIN_DLL / DLL and the -DTRACE, -DNFREQ... definitions come from the parent
# directories, so the EXPORT-marked ins_* functions are dllexport-ed here and
# dllimport-ed by the apps that link rtklib.

set(T_INS_SRCS
  insmath.c
  sins.c
  eskf.c
  imuread.c
  inscfg.c
  ins.c
)

add_library(t_ins OBJECT ${T_INS_SRCS})
set_property(TARGET t_ins PROPERTY C_STANDARD 99)
# built into the rtklib DLL, so the objects may not be position independent
set_property(TARGET t_ins PROPERTY POSITION_INDEPENDENT_CODE OFF)

target_include_directories(t_ins PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/src
)

if(MSVC)
  target_compile_definitions(t_ins PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

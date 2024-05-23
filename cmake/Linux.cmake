if (NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
  execute_process(COMMAND uname -s OUTPUT_VARIABLE OS_UNAME_SYSTEM)
  string(STRIP ${OS_UNAME_SYSTEM} OS_UNAME_SYSTEM)

  if (${OS_UNAME_SYSTEM} STREQUAL "Linux")
    add_definitions(-D_GNU_SOURCE)
    add_definitions(-DXOPEN_SOURCE=700)
    set(CMAKE_CXX_VISIBILITY_PRESET default)
  endif()
endif()


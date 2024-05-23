if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
  option (ENABLE_EMIT_LLVM "Emit object files as LLVM IR (with clang)" ON)
else()
  option (ENABLE_EMIT_LLVM "Emit object files as LLVM IR (with clang)" OFF)
endif()


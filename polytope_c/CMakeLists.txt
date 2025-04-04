# Build a configuration header file from our options.

configure_file(
  "${PROJECT_SOURCE_DIR}/polytope_c/polytope_c.h.in"
  "${PROJECT_BINARY_DIR}/polytope_c.h"
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(polytope_c polytope_c.cc polytope_plc.cc polytope_tessellator.cc
            polytope_tessellation.cc polytope_read_silo.cc polytope_write_silo.cc)
target_link_libraries(polytope_c polytope)

if (TESTING)
  add_subdirectory(tests)
endif()

# Library install targets
install (TARGETS polytope_c EXPORT polytope-targets DESTINATION lib)

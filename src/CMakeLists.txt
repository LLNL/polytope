# Include Silo stuff if we've found it.
if (SILO_FOUND)
  set(IO_SOURCES SiloWriter_2d.cc SiloWriter_3d.cc SiloReader_2d.cc SiloReader_3d.cc SiloUtils.cc)
endif()

# The tessellator sources we always build.
set(TESSELLATOR_SOURCES KeyTraits.cc predicates.cc PLC_CSG.cc QuantizedTessellation2d.cc QuantizedTessellation3d.cc
                        clipQuantizedTessellation2d.cc clipQuantizedTessellation3d.cc
                        snapToBoundary.cc)
set(CWD ${CMAKE_CURRENT_SOURCE_DIR})

# Check if we can include each of the serial tessellators:
if (USE_BOOST)
  list(APPEND TESSELLATOR_SOURCES MeshEditor.cc)
endif()

# TriangleTessellator and triangle source
if (HAVE_TRIANGLE)
  list(APPEND TESSELLATOR_SOURCES TriangleTessellator.cc)
  add_library(triangle STATIC  predicates.cc triangle.c)
  #set(TRIANGLE_LIB triangle PARENT_SCOPE)
  set(TRIANGLE_LIB triangle)
  set_target_properties(triangle PROPERTIES
    COMPILE_FLAGS "${COMPILE_FLAGS} -DTRILIBRARY -DCDT_ONLY -DANSI_DECLARATIONS -fPIC -DNO_TIMER"
    LINKER_LANGUAGE C)
endif()

# BoostTessellator
if (HAVE_BOOST_VORONOI)
  list(APPEND TESSELLATOR_SOURCES BoostTessellator.cc)
endif ()

# TetgenTessellator and tetgen source
if (HAVE_TETGEN)
  list(APPEND TESSELLATOR_SOURCES TetgenTessellator.cc)
  add_library(tetgen tetgen.cxx predicates.cc)
  #set(TETGEN_LIB tetgen PARENT_SCOPE)
  set(TETGEN_LIB tetgen)
  set_target_properties(${TETGEN_LIB} PROPERTIES
    COMPILE_FLAGS "${COMPILE_FLAGS} -DTETLIBRARY -Wuninitialized")
endif ()

# Check if we can incude the distributed tessellators:
if (HAVE_MPI)
  list(APPEND TESSELLATOR_SOURCES DistributedTessellator.cc SerialDistributedTessellator.cc)
endif()

# Add header files for IDE projects
file(GLOB HEADERS *.hh)
list(APPEND HEADERS "${PROJECT_BINARY_DIR}/polytope.hh")

add_library(polytopeC
  ${TESSELLATOR_SOURCES} ${IO_SOURCES} ErrorHandler.cc
  polytope_internal_abort.cc)

target_link_libraries(polytopeC
  ${TRIANGLE_LIB} ${TETGEN_LIB} ${SILO_LIBRARIES}
  ${HDF5_LIBRARIES} ${MPI_C_LIBRARIES}
  ${ZLIB_LIBRARIES})
# We must set the polytope python target to be "polytope"
# so the module is named polytope. So our C++ target is called
# polytopeC. But we will export it as polytope
set_target_properties(polytopeC PROPERTIES OUTPUT_NAME "polytope")
set_target_properties(polytopeC PROPERTIES EXPORT_NAME "polytope")
# Library install targets
install(TARGETS polytopeC EXPORT polytope-targets DESTINATION lib)
if (HAVE_TRIANGLE)
  install(TARGETS ${TRIANGLE_LIB} DESTINATION lib)
endif()
if (HAVE_TETGEN)
  install(TARGETS ${TETGEN_LIB} DESTINATION lib)
endif()

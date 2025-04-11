# ----------------------------------------------------------------------------------
# CMake file for Polytope's Python bindings using PYB11Generator
#   -Builds and compiles the Python bindings according to UsePyBindGen cmake module
# ----------------------------------------------------------------------------------

# The serial tessellators to be wrapped
set(SERIAL_TESSELLATORS "")
if(HAVE_TRIANGLE)
  set(SERIAL_TESSELLATORS ${SERIAL_TESSELLATORS} "TriangleTessellator2d")
endif()
if(HAVE_BOOST_VORONOI)
  set(SERIAL_TESSELLATORS ${SERIAL_TESSELLATORS} "BoostTessellator2d")
endif()
if(HAVE_TETGEN)
  set(SERIAL_TESSELLATORS ${SERIAL_TESSELLATORS} "TetgenTessellator3d")
endif()

# Generate the Tessellators.py script with the tessellators we're building
configure_file(Tessellators.py.in ${PROJECT_BINARY_DIR}/src/PYB11/Tessellators.py)

# libtriangle location
if(HAVE_TRIANGLE)
  set(TRIANGLE_STATIC_LIB "${PROJECT_BINARY_DIR}/src/libtriangle${CMAKE_STATIC_LIBRARY_SUFFIX}")
else()
  set(TRIANGLE_STATIC_LIB)
endif()

# libtetgen location
if(HAVE_TETGEN)
  set(TETGEN_STATIC_LIB "${PROJECT_BINARY_DIR}/src/libtetgen${CMAKE_STATIC_LIBRARY_SUFFIX}")
else()
  set(TETGEN_STATIC_LIB)
endif()

## Additional HDF5 link libraries
## NOTE:  This was added in to catch a few known build cases of static
##        silo and hdf5 on x86_64. Specifying the HDF5_ROOT will cause
##        Polytope to link against all the hdf5 libraries it finds.
## FIXME: A better solution is to allow users to specify a list of
##        additional link args when configuring Polytope.

# On some Linux machines we get messed up with the HDF5 libs including
# things like libdl and libm, so we try to cut down the list here to the
# first couple of elements.
list(LENGTH HDF5_LIBRARIES NUM_HDF5_LIBRARIES)
if (NUM_HDF5_LIBRARIES GREATER 2)
  list(GET HDF5_LIBRARIES 0 1 HDF5_LIB)
endif()

# Add our PYB11 module
PYB11Generator_add_module(polytope
  MODULE   polytope
  SOURCE   polytope_PYB11.py
  INCLUDES ${CMAKE_INSTALL_PREFIX}/include
  LINKS    polytopeC ${TRIANGLE_STATIC_LIB} ${TETGEN_STATIC_LIB} ${ZLIB_LIBRARIES} ${SILO_LIBRARIES} ${HDF5_LIBRARIES}
  INSTALL  ${POLYTOPE_SITE_PACKAGES_PATH}
  DEPENDS  polytopeC
  )

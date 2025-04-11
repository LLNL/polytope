# Minimum CMake version.
cmake_minimum_required (VERSION 3.1.0)

# Use XXX_ROOT variables to find third-party libs.
#cmake_policy(SET CMP0074 NEW)

# Set Polytope's module path.
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(Polytope)

# We require at least C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

##set(CMAKE_SKIP_RPATH ON)

# Project version numbers.
set(POLYTOPE_VERSION_MAJOR 1)
set(POLYTOPE_VERSION_MINOR 0)

# Set compilers. This must be done before enabling languages.
if (CC)
  set(CMAKE_C_COMPILER "${CC}")
  message(STATUS "C compiler is ${CMAKE_C_COMPILER}")
endif()
if (CXX)
  set(CMAKE_CXX_COMPILER "${CXX}")
  message(STATUS "C++ compiler is ${CMAKE_CXX_COMPILER}")
endif()

# Build everything as static libs.
option(BUILD_SHARED_LIBS "Build Polytope as shared library" OFF)

# Figure out the system type.
if (APPLE)
  set(SYS_FLAGS "-DAPPLE=1")
else ()
  if (LINUX)
  set(SYS_FLAGS "-DLINUX=1")
  endif ()
endif ()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SYS_FLAGS}")

# General compiler flags.
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare -ansi -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation")
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "XL")
  set(USING_XLC ON)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qPIC")
else()
  if (NOT MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation")
endif()
if (MSVC)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES)
endif (MSVC)

# Figure out MPI.
option(USE_MPI "Compile MPI parallel tessellator" OFF)
if (USE_MPI)
  message(FATAL ERROR "USE_MPI option is currently broken")
  find_package(MPI REQUIRED)
  message(STATUS "Using MPI for parallel tessellation.")
  set(HAVE_MPI ON)
  if(NOT MPI_CXX_INCLUDE_PATH STREQUAL "")
    include_directories("${MPI_CXX_INCLUDE_PATH}")
    find_program(SLURM_SRUN_COMMAND srun
      DOC "Path to Slurm's srun executable")
    if(SLURM_SRUN_COMMAND)
      set(MPIEXEC ${SLURM_SRUN_COMMAND})
      set(MPIEXEC_NUMPROC_FLAG -n)
      set(HAVE_MPIEXEC true)
      message("-- Slurm srun detected. Using ${MPIEXEC}")
    else()
      find_program(MPIRUN_COMMAND mpirun)
      if(MPIRUN_COMMAND)
        set(MPIEXEC mpirun)
        set(MPIEXEC_NUMPROC_FLAG -np)
        set(HAVE_MPIEXEC true)
      else()
        set(MPIEXEC "")
        set(MPIEXEC_NUMPROC_FLAG "")
        set(HAVE_MPIEXEC false)
      endif()
    endif()
  endif()
else()
  set(HAVE_MPI OFF)
endif()

# Find Boost.
option(USE_BOOST "Use Boost Voronoi tessellator" ON)

if (USE_BOOST)
  #set(Boost_ADDITIONAL_VERSIONS "1.51" "1.51.0" "1.52" "1.52.0" "1.53" "1.53.0" "1.60" )
  find_package(Boost 1.50 REQUIRED)
  if (Boost_FOUND)
    include_directories("${Boost_INCLUDE_DIR}")
    set(HAVE_BOOST ON)
    message(STATUS "Found Boost installation at ${Boost_INCLUDE_DIR}")
    if (Boost_MINOR_VERSION GREATER 51 AND
	EXISTS ${Boost_INCLUDE_DIR}/boost/polygon/voronoi.hpp)
      set(HAVE_BOOST_VORONOI ON)
      message(STATUS "Found Boost Voronoi in ${Boost_INCLUDE_DIR}/boost/polygon")
    else ()
      set(HAVE_BOOST_VORONOI OFF)
      message(STATUS "Boost Voronoi not found. BoostTessellator is disabled.")
    endif ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBoost_FOUND=1")
  else()
    set(HAVE_BOOST OFF)
  endif ()
else()
  set(HAVE_BOOST OFF)
endif ()

# Find Python
option(USE_PYTHON "Generate Python interface" OFF)
if (USE_PYTHON)
  message(STATUS "Looking for Python...")
  # If PYTHON_EXE is set through the bootstrap process, use it
  if (PYTHON_EXE)
    set(PYTHON_EXECUTABLE ${PYTHON_EXE})
    set(Python_EXECUTABLE ${PYTHON_EXE})
    set(Python3_EXECUTABLE ${PYTHON_EXE})
  endif()
  find_package(Python3 COMPONENTS Interpreter Development)
  if (${Python3_FOUND})
    set(HAVE_PYTHON ON)
    set(PYTHON_EXE ${Python3_EXECUTABLE})  # This is what PYB11Generator uses
    set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE}) # This is what pybind11 uses
    if (NOT PYB11GENERATOR_ROOT_DIR)
      set(PYB11GENERATOR_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern/PYB11Generator" CACHE PATH "")
    endif()
    if (NOT PYBIND11_ROOT_DIR)
      set(PYBIND11_ROOT_DIR "${PYB11GENERATOR_ROOT_DIR}/extern/pybind11" CACHE PATH "")
    endif()
    include(${PYB11GENERATOR_ROOT_DIR}/cmake/PYB11Generator.cmake)
    message(STATUS "Found Python version ${Python3_VERSION}")
    message(STATUS "Using Python executable ${Python3_EXECUTABLE}")
    message(STATUS "Using PYB11Generator ${PYB11GENERATOR_ROOT_DIR}")
    set(POLYTOPE_SITE_PACKAGES_PATH "lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages/polytope" CACHE PATH "")
  else()
    message(FATAL ERROR
      "USE_PYTHON set but no Python found. Try setting -DPYTHON_EXE=<local/python/install/bin/python>")
  endif()
else()
  set(HAVE_PYTHON OFF)
endif()

#Check for HDF5 and Silo
option(USE_SILO "Use SILO reader/writer" OFF)
set(ZLIB_LIBRARIES )
if (USE_SILO)
  # Find HDF5.
  set(HDF5_USE_STATIC_LIBRARIES 1)
  find_package(HDF5 REQUIRED)
  if (HDF5_FOUND)
    set(HAVE_HDF5 ON)
    # Get rid of'debug' and 'optimized' detritis.
    list(REMOVE_ITEM HDF5_LIBRARIES debug optimized)
    list(REMOVE_ITEM HDF5_INCLUDE_DIRS debug optimized)

    # Get rid of duplicates
    # NOTE: Taken from Will Dicharry's FindHDF5.cmake module. This ensures
    #       duplicates are removed from the beginning of the list so as
    #       not to mess with the unresolved symbols search at linking.
    list(REVERSE HDF5_LIBRARIES)
    list(REMOVE_DUPLICATES HDF5_LIBRARIES)
    list(REVERSE HDF5_LIBRARIES)
    message(STATUS "Found HDF5 headers in ${HDF5_INCLUDE_DIRS}")
    message(STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")

    # Remove "." from the include directories if it's in there.
    list(REMOVE_ITEM HDF5_INCLUDE_DIRS ".")

    # Do we have parallel HDF5?
    if (USE_MPI)
      if (HDF5_IS_PARALLEL)
        set(HAVE_PARALLEL_HDF5 ON)
      endif ()
    endif ()
  endif()

  # Try to find Silo. For simplicity, we assume it lives with HDF5.
  find_package(Silo REQUIRED)
  if (SILO_FOUND)
    set(HAVE_SILO ON)
    message(STATUS "Found Silo headers in ${SILO_INCLUDE_DIRS}")
    message(STATUS "Silo libraries: ${SILO_LIBRARIES}")
    include_directories("${SILO_INCLUDE_DIRS}")
    include_directories("${HDF5_INCLUDE_DIRS}")
  endif()

  # Silo also requires libz
  if (SILO_FOUND)
    find_package(ZLIB)
    include_directories("${ZLIB_INCLUDE_DIRS}")
    set(ZLIB_LIBRARIES "z")
  endif()

else()
  set(HAVE_HDF5 OFF)
  set(HAVE_SILO OFF)
endif()

# Enable Triangle if triangle.h and triangle.c are found.
if (EXISTS ${PROJECT_SOURCE_DIR}/src/triangle.h AND
    EXISTS ${PROJECT_SOURCE_DIR}/src/triangle.c AND
    Boost_FOUND)
  set(HAVE_TRIANGLE ON)
  message(STATUS "Triangle has been found! TriangleTessellator is enabled.")
else()
  set(HAVE_TRIANGLE OFF)
  message(STATUS "Triangle not found! TriangleTessellator is disabled.")
endif ()

# Enable Tetgen if the proper files are found.
if (EXISTS ${PROJECT_SOURCE_DIR}/src/tetgen.h AND
    EXISTS ${PROJECT_SOURCE_DIR}/src/tetgen.cxx)
  set(HAVE_TETGEN ON)
  message(STATUS "Tetgen has been found! TetgenTessellator is enabled.")
else()
  set(HAVE_TETGEN OFF)
  message(STATUS "Tetgen not found! TetgenTessellator is disabled.")
endif ()

# Build a configuration header file from our options.
configure_file(
  "${PROJECT_SOURCE_DIR}/src/polytope.hh.in"
  "${PROJECT_BINARY_DIR}/polytope.hh"
)

option(HEADER_ONLY "Header only installation" OFF)
if (NOT HEADER_ONLY)
  # Include the binary directory in the header file search path.
  include_directories(${PROJECT_BINARY_DIR})
  include_directories(${PROJECT_BINARY_DIR}/src)

  # Libraries

  # # Voro 2d/3d libraries.
  #include_directories(voro_2d)
  #add_subdirectory(voro_2d)
  #include_directories(voro_3d)
  #add_subdirectory(voro_3d)

  # Polytope proper.
  include_directories(src)
  add_subdirectory(src)

  # Tests directory
  option(TESTING "Generate tests" ON)
  if (TESTING)
    enable_testing()
    add_subdirectory(tests)
  endif ()

  # Polytope C library.
  if (BUILD_C_INTERFACE)
    message(STATUS "C interface (polytope_c) is enabled.")
    add_subdirectory(polytope_c)
  endif()

  # Python bindings.
  message(STATUS "Checking for python... ${HAVE_PYTHON}")
  if(HAVE_PYTHON)
    include_directories(src/PYB11)
    add_subdirectory(src/PYB11)
    message(STATUS "Adding src/PYB11")
  endif()
endif()

# Library install targets
install (FILES ${PROJECT_BINARY_DIR}/polytope.hh src/Tessellator.hh src/TessellatorInline.hh
               src/Tessellation.hh src/MeshEditor.hh src/ErrorHandler.hh
               src/PLC.hh src/PLC_CSG_2d.hh src/PLC_CSG_3d.hh
               src/ReducedPLC.hh src/simplifyPLCfacets.hh
               src/polytope_internal.hh src/QuantTessellation.hh
      	       src/QuantizedCoordinates.hh src/VoroPP_2d.hh src/VoroPP_3d.hh
               src/TriangleTessellator.hh src/BoostTessellator.hh
       	       src/BoostTessellatorTraits.hh src/TetgenTessellator.hh
      	       src/MeshEditor.hh src/polytope_internal.hh
               src/polytope_geometric_utilities.hh
               src/polytope_parallel_utilities.hh
               src/polytope_serialize.hh
      	       src/polytope_tessellator_utilities.hh
               src/polytope_write_OOGL.hh
      	       src/OrphanageBase.hh src/BoostOrphanage.hh src/Clipper2d.hh
               src/ErrorHandler.hh src/convexHull_2d.hh src/convexHull_3d.hh
               src/Point.hh src/KeyTraits.hh src/DimensionTraits.hh
               src/polytope_serialize.hh src/polytope_geometric_utilities.hh
               src/polytope_parallel_utilities.hh
      	       src/polytope_tessellator_utilities.hh
      	       src/PLC_CSG_2d.hh src/PLC_CSG_3d.hh src/PLC_Boost_2d.hh
               src/SiloWriter.hh src/SiloReader.hh src/polytope_write_OOGL.hh
               src/QuantizedTessellation2d.hh src/QuantizedTessellation3d.hh
               src/IntPointMap.hh src/clipQuantizedTessellation.hh
               src/removeElements.hh src/findBoundaryElements.hh
               src/snapToBoundary.hh src/makeBoxPLC.hh
         DESTINATION include/polytope)

# If we're parallel we have a few extra install items.
if (USE_MPI)
  install(FILES src/DistributedTessellator.hh
                src/SerialDistributedTessellator.hh
                src/checkDistributedTessellation.hh
          DESTINATION include/polytope)
endif ()

# Export package to allow it to be installed with find_package
install(EXPORT polytope-targets DESTINATION lib/cmake)
include(CMakePackageConfigHelpers)
set(POLYTOPE_INSTALL_INCLUDE_DIR "include" CACHE STRING "")
set(POLYTOPE_INSTALL_LIB_DIR "lib" CACHE STRING "")
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/polytope-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/polytope-config.cmake
  INSTALL_DESTINATION "lib"
  PATH_VARS
  POLYTOPE_INSTALL_INCLUDE_DIR
  POLYTOPE_INSTALL_LIB_DIR)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/polytope-config.cmake
  DESTINATION "lib/cmake")
# build a CPack driven installer package
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE
  "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${POLYTOPE_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${POLYTOPE_VERSION_MINOR}")
include(CPack)

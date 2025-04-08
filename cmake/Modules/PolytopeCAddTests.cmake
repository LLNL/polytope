# CMake macros for adding C tests into Polytope
# Consult test/CMakeLists.txt for usage
###############################################################

#--------------------------------------------------------------
# POLYTOPE_C_ADD_TEST
# Add a serial test to Polytope
#--------------------------------------------------------------
macro(POLYTOPE_C_ADD_TEST name dependency_list)
  # Every test links to the polytope_c and polytope libraries
  set(TEST_LINK_LIBRARIES polytope_c polytopeC)

  # Determine if you have the necessary components for the test
  set(BUILD_TEST true)
  foreach(_dependency ${dependency_list})
    set(DEP_NAME "HAVE_${_dependency}")
    if(NOT ${DEP_NAME})
      set(BUILD_TEST false)
    endif()
    # If using Tetgen, remember to link to its library
    if(${_dependency} EQUAL "TETGEN")
      set(APPEND EXTRA_LINK_LIBRARIES ${TETGEN_LIB})
    endif()
  endforeach()  

  # Build the test executable
  if(BUILD_TEST)
    set(TEST_NAME "test_${name}")
    add_executable(${TEST_NAME} "${TEST_NAME}.c")
    target_link_libraries(${TEST_NAME} ${TEST_LINK_LIBRARIES})

    # Special CTest run instructions
    add_test(${TEST_NAME} ${TEST_NAME})
    set_tests_properties(${TEST_NAME} 
      PROPERTIES PASS_REGULAR_EXPRESSION PASS
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
  endif()
endmacro()

#--------------------------------------------------------------
# POLYTOPE_C_ADD_DISTRIBUTED_TEST
# Add a parallel with custom run instructions
#--------------------------------------------------------------
macro(POLYTOPE_C_ADD_DISTRIBUTED_TEST name dependency_list procs)

  # Check for MPI and determine if you have the necessary
  # components to build the test.
  set(BUILD_TEST true)
  if(${HAVE_MPI} AND ${HAVE_MPIEXEC})
    # Every test links to the polytope library
    set(TEST_LINK_LIBRARIES polytopeC)
    foreach(_dependency ${dependency_list})
      set(DEP_NAME "HAVE_${_dependency}")
      if(NOT ${DEP_NAME})
	set(BUILD_TEST false)
      endif()
      # If using Tetgen, remember to link to its library
      if(${_dependency} EQUAL "TETGEN")
	set(APPEND EXTRA_LINK_LIBRARIES ${TETGEN_LIB})
      endif()
    endforeach()
  else()
    set(BUILD_TEST false)
  endif()

  # Build the test executable
  if(BUILD_TEST)
    set(TEST_NAME "test_${name}")
    add_executable(${TEST_NAME} "${TEST_NAME}.c")
    target_link_libraries(${TEST_NAME} ${TEST_LINK_LIBRARIES})
    # Special CTest run instructions
    #add_test(${TEST_NAME} ${TEST_NAME})
    foreach(proc ${procs})
      add_test(${TEST_NAME}_${proc}_proc 
	${MPIEXEC} 
	${MPIEXEC_NUMPROC_FLAG} 
	${proc}
	${MPIEXEC_PREFLAGS}
	${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}
	${MPIEXEC_POSTFLAGS})
      set_tests_properties(${TEST_NAME}_${proc}_proc 
	PROPERTIES PASS_REGULAR_EXPRESSION PASS
	WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
    endforeach()
  endif()
endmacro()

add_library(voro_3d cell.cc common.cc container.cc unitcell.cc v_compute.cc c_loops.cc
            v_base.cc wall.cc pre_container.cc container_prd.cc)

# Library install targets
install (TARGETS voro_3d DESTINATION lib EXPORT polytope-targets)

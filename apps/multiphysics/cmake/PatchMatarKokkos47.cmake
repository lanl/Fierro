# A DualView whose memory space is host-accessible must share its allocation.
# create_mirror_view_and_copy aliases host memory, and copies device memory.
set(_header "${MATAR_SOURCE_DIR}/src/include/kokkos_types.h")
file(READ "${_header}" _content)
set(_old "mystrides_host_ = typename Strides1D::t_host(\"host_strides\", dims_[0]);
    //requires host synchronization before building dual view wrapper
    Kokkos::deep_copy(mystrides_host_, strides_array.get_kokkos_view());")
set(_new "mystrides_host_ = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace{}, strides_array.get_kokkos_view());")
string(FIND "${_content}" "${_old}" _old_pos)
string(FIND "${_content}" "${_new}" _new_pos)
if(_old_pos GREATER_EQUAL 0)
  string(REPLACE "${_old}" "${_new}" _content "${_content}")
  file(WRITE "${_header}" "${_content}")
elseif(_new_pos LESS 0)
  message(FATAL_ERROR "Unrecognized MATAR DualView constructors in ${_header}")
endif()

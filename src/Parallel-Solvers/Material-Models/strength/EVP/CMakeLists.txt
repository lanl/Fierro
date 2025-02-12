set(SRC_Files EVPStrengthModel.cpp chg_basis.cpp evpal.cpp user_mat.cpp)
add_library(EVPStrengthModel OBJECT ${SRC_Files})
target_include_directories(EVPStrengthModel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(EVPStrengthModel Elements)
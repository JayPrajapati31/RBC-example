# Sriramajayam

cmake_minimum_required(VERSION 3.12)

project(rbc-axisymm-examples
  VERSION 1.0
  DESCRIPTION "Examples using rbc with nlopt")

foreach(TARGET
    radial_RHBC
    radial_RHBC_1
    radial_RHBC_2
    radial_RHBC_3)

# find dependencies
find_package(rbc REQUIRED)

list(APPEND COMPILE_FEATURES
  "cxx_constexpr"
  "cxx_auto_type"
  "cxx_decltype"
  "cxx_decltype_auto"
  "cxx_enum_forward_declarations"
  "cxx_inline_namespaces"
  "cxx_lambdas"
  "cxx_nullptr"
  "cxx_override"
  "cxx_range_for")

  #set(TARGET "multi_pass")

# add executable
add_executable(${TARGET} ${TARGET}.cpp customFunc.cpp ReadJSON.cpp)

# include headers
target_include_directories(${TARGET} PUBLIC ${RBC_INCLUDE_DIRS} $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# link directories
target_link_directories(${TARGET} PUBLIC ${RBC_LIBRARY_DIRS})
target_link_libraries(${TARGET} PUBLIC ${RBC_LIBRARIES} ${CMAKE_DL_LIBS} NLopt::nlopt)

# choose appropriate flags
target_compile_features(${TARGET} PUBLIC ${COMPILE_FEATURES})

endforeach()

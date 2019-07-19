set(_common_libraries
  ActsCore
  ACTFramework
  ACTFWExamplesCommon
  ActsFrameworkGenerators
  ActsFrameworkPythia8
  ActsFrameworkTruthTracking
  ACTFWCsvPlugin
  ACTFWRootPlugin
  ${Boost_LIBRARIES})

# Fitter example
file(GLOB_RECURSE src_fit_files "src/VertexFit*.cpp")
add_executable(ACTFWVertexFitterExample ${src_fit_files})
target_include_directories(ACTFWVertexFitterExample PRIVATE
  ${Boost_INCLUDE_DIRS})
target_link_libraries(ACTFWVertexFitterExample PRIVATE
  ${_common_libraries})
install(
  TARGETS ACTFWVertexFitterExample
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})


# Finder example
file(GLOB_RECURSE src_find_files "src/VertexFind*.cpp")
add_executable(ACTFWVertexFinderExample ${src_find_files})
target_include_directories(ACTFWVertexFinderExample PRIVATE
  ${Boost_INCLUDE_DIRS})
target_link_libraries(ACTFWVertexFinderExample PRIVATE
  ${_common_libraries})
install(
  TARGETS ACTFWVertexFinderExample
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

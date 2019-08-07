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

  # Writer example
file(GLOB_RECURSE src_wr_files "src/VertexAndTrackWriterExample.cpp")
add_executable(ACTFWVertexWriterExample ${src_wr_files})
target_include_directories(ACTFWVertexWriterExample PRIVATE
  ${Boost_INCLUDE_DIRS})
target_link_libraries(ACTFWVertexWriterExample PRIVATE
  ${_common_libraries})
install(
  TARGETS ACTFWVertexWriterExample
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

set(sources src/VertexAndTracksReaderAndFinderExample.cpp
            src/VertexFindingAlgorithm.cpp)

  # Reader example
add_executable(ACTFWVertexReaderExample ${sources})
target_include_directories(ACTFWVertexReaderExample PRIVATE
  ${Boost_INCLUDE_DIRS})
target_link_libraries(ACTFWVertexReaderExample PRIVATE
  ${_common_libraries})
install(
  TARGETS ACTFWVertexReaderExample
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

set(_common_libraries
  ActsCore
  ACTFramework
  ACTFWPropagation
  ACTFWExamplesCommon
  ACTFWBFieldPlugin
  ACTFWObjPlugin
  ACTFWJsonPlugin
  ACTFWRootPlugin
  ACTFWPlyPlugin
  ${Boost_LIBRARIES})

# Generic detector
add_executable(ACTFWPlyWriterExample src/PlyWriterExample.cpp)
target_include_directories(ACTFWPlyWriterExample PRIVATE
  ${Boost_INCLUDE_DIRS})
target_link_libraries(ACTFWPlyWriterExample PRIVATE
  ${_common_libraries}
  ACTFWGenericDetector)
install(
  TARGETS ACTFWPlyWriterExample
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
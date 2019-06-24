set(_sources
  src/Barcode/BarcodeSvc.cpp
  src/Framework/BareAlgorithm.cpp
  src/Framework/Sequencer.cpp
  src/Random/RandomNumbersSvc.cpp
  src/Utilities/Paths.cpp
  src/Utilities/Helpers.cpp
  src/Validation/EffPlotTool.cpp
  src/Validation/ResPlotTool.cpp)

add_library(ACTFramework SHARED ${_sources})
target_include_directories(
  ACTFramework
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> ${ROOT_INCLUDE_DIRS}
  PRIVATE ${TBB_INCLUDE_DIRS})
target_link_libraries(
  ACTFramework
  PUBLIC ActsCore ${ROOT_LIBRARIES}
  PRIVATE ${TBB_LIBRARIES} dfelibs)
install(TARGETS ACTFramework LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

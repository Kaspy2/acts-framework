add_library(ACTFramework SHARED
  src/Framework/BareAlgorithm.cpp
  src/Framework/BareService.cpp
  src/Framework/RandomNumbers.cpp
  src/Framework/Sequencer.cpp
  src/Utilities/Paths.cpp)
target_include_directories(
  ACTFramework
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> ${ROOT_INCLUDE_DIRS}
  PRIVATE ${TBB_INCLUDE_DIRS})
target_link_libraries(
  ACTFramework
  PUBLIC ActsCore ${ROOT_LIBRARIES}
  PRIVATE ${TBB_LIBRARIES} dfelibs)
install(TARGETS ACTFramework LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

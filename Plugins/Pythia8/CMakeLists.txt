add_library(ActsFrameworkPythia8 SHARED
  ACTFW/Generators/Pythia8ProcessGenerator.cpp
  ACTFW/Options/Pythia8Options.cpp)

target_include_directories(ActsFrameworkPythia8 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_include_directories(ActsFrameworkPythia8 PUBLIC
  ${Boost_INCLUDE_DIRS})
target_link_libraries(ActsFrameworkPythia8 PUBLIC
  ActsCore
  ActsFrameworkGenerators
  ACTFramework
  Pythia8
  ${Boost_LIBRARIES})

install(TARGETS ActsFrameworkPythia8 LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY ACTFW DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
set(Generators_SOURCES
  ACTFW/Generators/EventGenerator.cpp
  ACTFW/Generators/ParametricProcessGenerator.cpp
  ACTFW/Generators/ParticleSelector.cpp
  ACTFW/Options/ParticleGunOptions.cpp)

add_library(ActsFrameworkGenerators SHARED ${Generators_SOURCES})
target_include_directories(ActsFrameworkGenerators
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
target_link_libraries(ActsFrameworkGenerators
  PUBLIC ActsCore ACTFramework ${Boost_PROGRAM_OPTIONS_LIBRARY})

install(TARGETS ActsFrameworkGenerators LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})


set(_common_libraries
  ActsCore
  ACTFramework
  ACTFWExamplesCommon
  ${Boost_LIBRARIES})

set(_plugin_libraries
  ACTFWRootPlugin
  ACTFWJsonPlugin
  ACTFWBFieldPlugin)

add_executable(ACTFWGenericMaterialValidationExample src/GenericMaterialValidation.cpp)
target_include_directories(ACTFWGenericMaterialValidationExample PRIVATE 
  ${Boost_INCLUDE_DIRS})
target_link_libraries(ACTFWGenericMaterialValidationExample PRIVATE 
  ${_common_libraries} ${_plugin_libraries})
target_link_libraries(ACTFWGenericMaterialValidationExample PRIVATE 
  ACTFWPropagation ACTFWGenericDetector)
install(TARGETS ACTFWGenericMaterialValidationExample RUNTIME DESTINATION 
  ${CMAKE_INSTALL_BINDIR})

add_executable(ACTFWGenericMaterialMappingExample src/GenericMaterialMapping.cpp)
target_include_directories(ACTFWGenericMaterialMappingExample PRIVATE 
  ${Boost_INCLUDE_DIRS})
target_link_libraries(ACTFWGenericMaterialMappingExample PRIVATE 
  ${_common_libraries} ${_plugin_libraries})
target_link_libraries(ACTFWGenericMaterialMappingExample PRIVATE 
  ACTFWMaterialMapping ACTFWGenericDetector)
install(TARGETS ACTFWGenericMaterialMappingExample RUNTIME DESTINATION 
  ${CMAKE_INSTALL_BINDIR})


if (USE_GEANT4 AND USE_DD4HEP)
  include(${Geant4_USE_FILE})
  
  add_executable(ACTFWGeantinoRecordingExample src/GeantinoRecordingExample.cpp)
  target_include_directories(ACTFWGeantinoRecordingExample PRIVATE 
    ${Boost_INCLUDE_DIRS} ${DD4hep_INCLUDE_DIRS})
  target_link_libraries(ACTFWGeantinoRecordingExample PRIVATE 
   ACTFWExamplesCommon ACTFWMaterialMapping ACTFWDD4hepG4Plugin ACTFWRootPlugin)
  target_link_libraries(ACTFWGeantinoRecordingExample PRIVATE 
   ACTFWDD4hepDetector ${DD4hep_LIBRARIES} ${Geant4_LIBRARIES} ${Boost_LIBRARIES})
   
  install(TARGETS ACTFWGeantinoRecordingExample RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  
  add_executable(ACTFWParticleRecordingExample src/ParticleRecordingExample.cpp)
  target_include_directories(ACTFWParticleRecordingExample PRIVATE 
    ${Boost_INCLUDE_DIRS} ${DD4hep_INCLUDE_DIRS})
  target_link_libraries(ACTFWParticleRecordingExample PRIVATE 
   ACTFWExamplesCommon ACTFWMaterialMapping ACTFWDD4hepG4Plugin)
  target_link_libraries(ACTFWParticleRecordingExample PRIVATE 
   ACTFWDD4hepDetector ${DD4hep_LIBRARIES} ${Geant4_LIBRARIES} ${Boost_LIBRARIES})
   
  install(TARGETS ACTFWParticleRecordingExample RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

set(_common_libraries
  ActsCore
  FatrasCore
  ACTFramework
  ACTFWBFieldPlugin
  ACTFWCsvPlugin
  ACTFWDigitization
  ACTFWExamplesCommon
  ACTFWFatras
  ACTFWFitting
  ACTFWObjPlugin
  ACTFWRootPlugin
  ActsFrameworkGenerators
  ActsFrameworkPythia8)

# Generic detector
add_executable(ACTFWGenericFatrasExample src/GenericFatrasExample.cpp)
target_include_directories(ACTFWGenericFatrasExample PRIVATE
  ${Boost_INCLUDE_DIRS})
target_link_libraries(ACTFWGenericFatrasExample PRIVATE
  ${_common_libraries}
  ACTFWGenericDetector)
install(
  TARGETS ACTFWGenericFatrasExample
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Generic detector with IOV based alignment
add_executable(ACTFWAlignedFatrasExample src/AlignedFatrasExample.cpp)
target_include_directories(ACTFWAlignedFatrasExample PRIVATE
  ${Boost_INCLUDE_DIRS})
target_link_libraries(ACTFWAlignedFatrasExample PRIVATE
  ${_common_libraries}
  ACTFWContextualDetector)
install(
  TARGETS ACTFWAlignedFatrasExample
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Generic detector with Payload based alignment
add_executable(ACTFWPayloadFatrasExample src/PayloadFatrasExample.cpp)
target_include_directories(ACTFWPayloadFatrasExample PRIVATE
  ${Boost_INCLUDE_DIRS})
target_link_libraries(ACTFWPayloadFatrasExample PRIVATE
  ${_common_libraries}
  ACTFWContextualDetector)
install(
  TARGETS ACTFWPayloadFatrasExample
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# TGEO based detector
if (USE_TGEO)
  add_executable(ACTFWTGeoFatrasExample src/TGeoFatrasExample.cpp)
  target_include_directories(ACTFWTGeoFatrasExample PRIVATE
    ${Boost_INCLUDE_DIRS})
  target_link_libraries(ACTFWTGeoFatrasExample PRIVATE
    ${_common_libraries}
    ACTFWTGeoDetector)
  install(
    TARGETS ACTFWTGeoFatrasExample
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# DD4hep detector
if (USE_DD4HEP)
  add_executable(ACTFWDD4hepFatrasExample src/DD4hepFatrasExample.cpp)
  target_include_directories(ACTFWDD4hepFatrasExample PRIVATE
    ${Boost_INCLUDE_DIRS}
    ${DD4hep_INCLUDE_DIRS})
  target_link_libraries(ACTFWDD4hepFatrasExample PRIVATE
    ${_common_libraries}
    ACTFWDD4hepDetector
    ${DD4hep_LIBRARIES})
  install(
    TARGETS ACTFWDD4hepFatrasExample
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

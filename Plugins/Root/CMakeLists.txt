set(sources src/RootMaterialDecorator.cpp
            src/RootMaterialWriter.cpp            
            src/RootMaterialTrackReader.cpp
            src/RootMaterialTrackWriter.cpp
            src/RootPlanarClusterWriter.cpp
            src/RootParticleWriter.cpp
            src/RootPropagationStepsWriter.cpp
            src/RootSimHitWriter.cpp
            src/RootTrackWriter.cpp
            src/RootPerformanceValidation.cpp
          )

add_library(ACTFWRootPlugin SHARED ${sources})
target_include_directories(
  ACTFWRootPlugin
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> ${ROOT_INCLUDE_DIRS})
target_link_libraries(
  ACTFWRootPlugin
  PUBLIC ActsCore ActsDigitizationPlugin IdentificationPlugin
         ACTFramework ACTFWPropagation ${ROOT_LIBRARIES} Threads::Threads)
install(TARGETS ACTFWRootPlugin LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

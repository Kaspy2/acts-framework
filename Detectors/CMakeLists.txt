add_subdirectory(ContextualDetector)
add_subdirectory_if(DD4hepDetector USE_DD4HEP)
add_subdirectory(GenericDetector)
add_subdirectory_if(TGeoDetector USE_TGEO)

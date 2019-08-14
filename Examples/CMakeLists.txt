# shared code
add_subdirectory(Common)

# tools
add_subdirectory(BField)
add_subdirectory(EventGenerator)
add_subdirectory_if(Fatras USE_PYTHIA8)
add_subdirectory(Geometry)
add_subdirectory(HelloWorld)
add_subdirectory_if(HepMCExample USE_HEPMC3)
add_subdirectory(MaterialMapping)
add_subdirectory(Propagation)
add_subdirectory_if(Vertexing USE_PYTHIA8)
add_subdirectory_if(HepMCExample USE_HEPMC3)

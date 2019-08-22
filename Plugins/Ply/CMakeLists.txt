file(GLOB_RECURSE src_files "src/*.cpp" "include/*.*pp")

add_library(ACTFWPlyPlugin SHARED ${src_files})
target_include_directories(ACTFWPlyPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(ACTFWPlyPlugin PUBLIC ActsCore ACTFramework)
target_link_libraries(ACTFWPlyPlugin PUBLIC Threads::Threads)

install(TARGETS ACTFWPlyPlugin LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ACTFW DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

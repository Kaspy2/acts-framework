add_executable(ACTFWHelloWorldExample
  HelloWorld.cpp
  HelloLoggerAlgorithm.cpp
  HelloRandomAlgorithm.cpp
  HelloService.cpp
  HelloWhiteBoardAlgorithm.cpp)
target_include_directories(ACTFWHelloWorldExample
  PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(ACTFWHelloWorldExample
  PRIVATE ACTFramework ACTFWExamplesCommon ActsCore ${Boost_LIBRARIES})
install(
  TARGETS ACTFWHelloWorldExample RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

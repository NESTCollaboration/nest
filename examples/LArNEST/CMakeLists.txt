# add_executable(bareLArNEST bareLArNEST.cpp)
# target_link_libraries(bareLArNEST PUBLIC NEST::Core)

add_executable(LegacyLArNESTBenchmarks LegacyLArNESTBenchmarks.cpp)
target_link_libraries(LegacyLArNESTBenchmarks PUBLIC NEST::Core)

add_executable(LArNESTBenchmarks LArNESTBenchmarks.cpp)
target_link_libraries(LArNESTBenchmarks PUBLIC NEST::Core)

add_executable(LArNESTNeutronCapture LArNESTNeutronCapture.cpp)
target_link_libraries(LArNESTNeutronCapture PUBLIC NEST::Core)
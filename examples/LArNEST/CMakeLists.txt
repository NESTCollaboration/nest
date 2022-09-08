# add_executable(bareLArNEST bareLArNEST.cpp)
# target_link_libraries(bareLArNEST PUBLIC NEST::Core)

add_executable(LegacyLArNESTBenchmarks LegacyLArNESTBenchmarks.cpp)
target_link_libraries(LegacyLArNESTBenchmarks PUBLIC NEST::Core)

add_executable(LArNESTMeanYieldsBenchmarks LArNESTMeanYieldsBenchmarks.cpp)
target_link_libraries(LArNESTMeanYieldsBenchmarks PUBLIC NEST::Core)

add_executable(LArNESTFluctuationBenchmarks LArNESTFluctuationBenchmarks.cpp)
target_link_libraries(LArNESTFluctuationBenchmarks PUBLIC NEST::Core)

add_executable(LArNESTNeutronCapture LArNESTNeutronCapture.cpp)
target_link_libraries(LArNESTNeutronCapture PUBLIC NEST::Core)
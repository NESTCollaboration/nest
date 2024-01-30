
set(MIN_ROOT_VERSION 6.12)

add_executable(GenerateGasTable GenerateGarfieldGasTableForLiquidNoble.cpp)
target_link_libraries(GenerateGasTable PUBLIC NESTCore)
set(GARFIELD GenerateGasTable)

install(
    TARGETS ${GARFIELD}        
    RUNTIME DESTINATION bin
)









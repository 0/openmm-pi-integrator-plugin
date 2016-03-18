set(target PiIntegratorReference)

file(GLOB source_files src/*.cpp)
file(GLOB include_files include/*.h)

add_library(${target} SHARED ${source_files})
include_directories(BEFORE "${CMAKE_SOURCE_DIR}/openmmapi/include")
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(${target} OpenMM)
install(TARGETS ${target} LIBRARY DESTINATION lib/plugins)

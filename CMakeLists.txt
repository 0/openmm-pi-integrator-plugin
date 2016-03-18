cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

set(OPENMM_DIR /usr/local/openmm CACHE PATH "Where OpenMM is installed")
set(CMAKE_INSTALL_PREFIX ${OPENMM_DIR} CACHE PATH "Where to install the plugin" FORCE)

include_directories("${OPENMM_DIR}/include")
link_directories("${OPENMM_DIR}/lib" "${OPENMM_DIR}/lib/plugins")

### Subdirectories.
add_subdirectory(platforms/reference)
add_subdirectory(python)

### Main library.
set(target PiIntegrator)

file(GLOB source_files openmmapi/src/*.cpp)
file(GLOB include_files openmmapi/include/*.h)

add_library(${target} SHARED ${source_files})
include_directories(BEFORE openmmapi/include)
target_link_libraries(${target} OpenMM)
install(TARGETS ${target} LIBRARY DESTINATION lib)
install(FILES ${include_files} DESTINATION include)

find_package(SWIG 3 REQUIRED)

set(interface_file PiIntegratorPlugin.i)
set(wrap_file PiIntegratorPluginWrapper.cpp)

add_custom_command(
    OUTPUT ${wrap_file}
    COMMAND ${SWIG_EXECUTABLE} "-I${OPENMM_DIR}/include" -c++ -python -o ${wrap_file} ${interface_file}
    DEPENDS ${interface_file}
    VERBATIM
    )

set(PY_PLUGIN_NAME PiIntegratorPlugin)
set(PY_PLUGIN_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/openmmapi/include")
set(PY_PLUGIN_LIBRARY_DIR ${CMAKE_BINARY_DIR})
configure_file(setup.template.py setup.py)

add_custom_target(
    PythonInstall
    COMMAND python setup.py build
    COMMAND python setup.py install
    DEPENDS ${wrap_file}
    VERBATIM
    )

#!/usr/bin/env python

from distutils.core import Extension, setup
from os.path import join
import platform

openmm_dir = '@OPENMM_DIR@'
plugin_include_dir = '@PY_PLUGIN_INCLUDE_DIR@'
plugin_library_dir = '@PY_PLUGIN_LIBRARY_DIR@'

extra_link_args = []

if platform.system() == 'Darwin':
    extra_link_args.append('-Wl,-rpath,' + join(openmm_dir, 'lib'))

module = Extension('_@PY_PLUGIN_NAME@',
                   sources=['@wrap_file@'],
                   libraries=['OpenMM', 'PiIntegrator'],
                   include_dirs=[join(openmm_dir, 'include'), plugin_include_dir],
                   library_dirs=[join(openmm_dir, 'lib'), plugin_library_dir],
                   extra_link_args=extra_link_args,
                   )

setup(name='@PY_PLUGIN_NAME@',
      version='0.1',
      author='Dmitri Iouchtchenko',
      description='Path integral integrator plugin for OpenMM.',
      ext_modules=[module],
      py_modules=['@PY_PLUGIN_NAME@'],
      )

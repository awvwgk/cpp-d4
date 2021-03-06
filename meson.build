# This file is part of cpp-d4.
#
# Copyright (C) 2019-2020 Sebastian Ehlert
#
# cpp-d4 is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cpp-d4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with cpp-d4.  If not, see <https://www.gnu.org/licenses/>.

project(
  'cpp-d4',
  'cpp',
  version: '1.0',
  license: 'LGPL3',
  meson_version: '>=0.45',
)

cxx = meson.get_compiler('cpp')

deps = []
foreach lib: ['cblas', 'lapacke']
  this_dep = dependency(lib)
  if not this_dep.found()
    this_dep = cxx.find_library(lib)
  endif
  deps += this_dep
endforeach

srcs = files(
  'src/dftd_dispersion.cpp',
  'src/dftd_eeq.cpp',
  'src/dftd_ncoord.cpp',
  'src/dftd_damping.cpp',
)

cpp_d4_inc = include_directories('include')

cpp_d4_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  include_directories: cpp_d4_inc,
  dependencies: deps,
  install: true,
)

install_headers(
  'include/dftd_dispersion.h',
  'include/dftd_eeq.h',
  'include/dftd_ncoord.h',
  'include/dftd_parameters.h',
  'include/dftd_damping.h',
  'include/dftd_matrix.h',
  'include/dftd_geometry.h',
)

cpp_d4_dep = [
  declare_dependency(
    link_with: cpp_d4_lib,
    include_directories: cpp_d4_inc,
  ),
  deps,
]

executable(
  meson.project_name(),
  sources: files(
    'src/program_dftd.cpp',
    'src/dftd_readxyz.cpp',
  ),
  dependencies: cpp_d4_dep,
  install: true,
)

/* This file is part of cpp-d4.
 *
 * Copyright (C) 2019 Sebastian Ehlert, Marvin Friede
 *
 * cpp-d4 is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * cpp-d4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with cpp-d4.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef UTIL_H
#define UTIL_H

#include <string>

#include <dftd_geometry.h>

extern int get_molecule(
  int n,
  const char atoms[][3],
  const double coord[],
  dftd4::TMolecule &mol
);

extern bool check(
  double actual,
  double expected,
  double epsilon = 1e-12,
  bool rel = false
);
extern bool
  check(float actual, float expected, float epsilon = 1e-6, bool rel = false);

extern void
  print_fail(const char specifier[32], double obtained, double expected);

extern int element(const std::string &sym);

#endif // UTIL_H

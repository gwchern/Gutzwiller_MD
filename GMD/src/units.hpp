//
//  units.hpp
//  GQMD-Hubbard-liquid
//
//  Created by Gia-Wei Chern on 1/12/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#ifndef units_h
#define units_h

// length
constexpr double angstrom = 1.0;
constexpr double bohr = 0.5292 * angstrom; // 0.5291772109217

// energy
constexpr double eV = 1.0;
constexpr double rydberg = 13.605804 * eV; // 13.6056925330
constexpr double hartree = 2 * rydberg;

// mass
constexpr double amu = 1.0;
constexpr double massSi = 28.0851 * amu;
constexpr double mass_e = amu / 1822.888;

// derived units
constexpr double hbar = sqrt(hartree * mass_e) * bohr; //0.0646572

// time, derived unit, because length energy and mass are specified
constexpr double femtosecond = /* sqrt(amu/eV) * angstrom */ 1. / 10.1805054836529;

// charge
constexpr double charge_e = 1.0;

// temperature
constexpr double kelvin = 1.0;
constexpr double kB = 8.6173323849609e-5 * eV / kelvin; //8.61733e-05,    kT = 0.1 ~ 1160.5 K

#endif /* units_h */

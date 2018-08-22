// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdlib>
#include <iostream>
#include <utility>
#include "ACTFW/ParticleGun/ParticleGun.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "Acts/Utilities/Units.hpp"

namespace po = boost::program_options;

namespace au = Acts::units;

namespace FW {

namespace Options {

  /// The particle gun options, the are prefixes with gp
  ///
  /// @tparam aopt_t Type of the options object
  ///
  /// @param opt the options object to be parsed
  template <typename aopt_t>
  void
  addParticleGunOptions(aopt_t& opt)
  {
    opt.add_options()("pg-nparticles",
                      po::value<size_t>()->default_value(100.),
                      "number of particles.")(
        "pg-pdg",
        po::value<int>()->default_value(13),
        "PDG number of the particle, will be adjusted for charge flip.")(
        "pg-mass",
        po::value<double>()->default_value(105.),
        "mass of the particle in [MeV]")(
        "pg-charge",
        po::value<double>()->default_value(-1.),
        "charge of the particle in [e]")(
        "pg-chargeflip",
        po::value<bool>()->default_value(true),
        "flip the charge (and change PDG accordingly).")(
        "pg-d0-range",
        po::value<read_range>()->multitoken()->default_value({0., 0.}),
        "range in which the d0 parameter is simulated in [mm]. Please hand"
        "over by simply seperating the values by space")(
        "pg-z0-range",
        po::value<read_range>()->multitoken()->default_value({0., 0.}),
        "range in which the z0 parameter is simulated in [mm]. Please hand"
        "over by simply seperating the values by space")(
        "pg-phi-range",
        po::value<read_range>()->multitoken()->default_value({-M_PI, M_PI}),
        "range in which the phi0 parameter is simulated. Please hand over by "
        "simply seperating the values by space")(
        "pg-eta-range",
        po::value<read_range>()->multitoken()->default_value({-4., 4.}),
        "range in which the eta parameter is simulated. Please hand over by "
        "simply seperating the values by space")(
        "pg-pt-range",
        po::value<read_range>()->multitoken()->default_value({0.1, 1e3}),
        "range in which the pt in [GeV] parameter is simulated. Please hand "
        "over by simply seperating the values by space");
  }

  /// Read the particle gun options and return a Config file
  /// 
  /// @tparam amap_t Type of the map object
  ///
  /// @param vm The map object
  template <typename amap_t>
  FW::ParticleGun::Config
  readParticleGunConfig(const amap_t& vm)
  {
    // read the reange as vector (missing istream for std::array)
    auto d0r  = vm["pg-d0-range"].template as<read_range>();
    auto z0r  = vm["pg-z0-range"].template as<read_range>();
    auto phir = vm["pg-phi-range"].template as<read_range>();
    auto etar = vm["pg-eta-range"].template as<read_range>();
    auto ptr  = vm["pg-pt-range"].template as<read_range>();
    // particle gun as generator
    FW::ParticleGun::Config particleGunConfig;
    particleGunConfig.nParticles = vm["pg-nparticles"].template as<size_t>();
    particleGunConfig.d0Range    = {{d0r[0] * au::_mm, d0r[1] * au::_mm}};
    particleGunConfig.z0Range    = {{z0r[0] * au::_mm, z0r[1] * au::_mm}};
    particleGunConfig.phiRange   = {{phir[0], phir[1]}};
    particleGunConfig.etaRange   = {{etar[0], etar[1]}};
    particleGunConfig.ptRange    = {{ptr[0] * au::_GeV, ptr[1] * au::_GeV}};
    particleGunConfig.mass   = vm["pg-mass"].template as<double>() * au::_MeV;
    particleGunConfig.charge = vm["pg-charge"].template as<double>() * au::_e;
    particleGunConfig.randomCharge = vm["pg-chargeflip"].template as<bool>();
    particleGunConfig.pID          = vm["pg-pdg"].template as<int>();
    // return the config object
    return particleGunConfig;
  }
}
}

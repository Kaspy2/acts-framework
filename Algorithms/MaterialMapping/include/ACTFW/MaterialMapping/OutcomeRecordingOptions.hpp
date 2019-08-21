// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdlib>
#include <iostream>
#include <utility>
#include "ACTFW/Utilities/Options.hpp"
#include "Acts/Utilities/Units.hpp"

namespace po = boost::program_options;

namespace au = Acts::units;

namespace FW {

namespace Options {

  /// the particle gun options, the are prefixes with gp
  template <typename aopt_t>
  void
  addOutcomeRecordingOptions(aopt_t& opt)
  {
    opt.add_options()(
        "ammo",
        po::value<std::string>()->default_value("proton"),
        "Ammo name of the particle gun")(
        "tracks-per-event",
        po::value<size_t>()->default_value(1),
        "Number of tracks produced per event")(
        "energy",
        po::value<double>()->default_value(1000.),
        "Initial energy of the primary particle in MeV.")(
        "lockangle",
        po::value<bool>()->default_value(false),
        "Boolean flag to fix a certain angular direction of the primary particle.")(
        "phi",
        po::value<double>()->default_value(0.),
        "Phi angle of the primary particle.")(
        "Theta",
        po::value<double>()->default_value(0.5 * M_PI),
        "Theta angle of the primary particle.")(
        "lockpos",
        po::value<bool>()->default_value(false),
        "Boolean flag to fix a certain vertex position of the primary particle.")(
        "x",
        po::value<double>()->default_value(0.),
        "x position of the primary vertex")(
        "y",
        po::value<double>()->default_value(0.),
        "y position of the primary vertex")(
        "z",
        po::value<double>()->default_value(0.),
        "z position of the primary vertex");
  }
  
  /// read the particle gun options and return a Config file
  template <typename amap_t>
  FW::OutcomeRecording::Config
  readOutcomeRecordingConfig(const amap_t& vm)
  {
    double x = 0., y = 0., z = 0.;
    
    FW::OutcomeRecording::Config config;
    
    config.particleName = vm["ammo"].template as<std::string>();
    config.tracksPerEvent = vm["tracks-per-event"].template as<size_t>();
	config.energy = vm["energy"].template as<double>();
	config.lockAngle = vm["lockangle"].template as<bool>();
	config.phi = vm["phi"].template as<double>();
	config.theta = vm["theta"].template as<double>();
	config.lockPosition = vm["lockpos"].template as<bool>();
	x = vm["x"].template as<double>();
	y = vm["y"].template as<double>();
	z = vm["z"].template as<double>();
	config.pos = {x, y, z};
    
    return config;
  }
}
}
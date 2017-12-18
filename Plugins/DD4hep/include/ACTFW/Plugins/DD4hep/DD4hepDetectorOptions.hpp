// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTFW_OPTIONS_DD4HEPDETECTOROPTIONS_HPP
#define ACTFW_OPTIONS_DD4HEPDETECTOROPTIONS_HPP

#include <cstdlib>
#include <iostream>
#include <utility>
#include "ACTFW/Plugins/DD4hep/GeometryService.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "Acts/Utilities/Units.hpp"

namespace po = boost::program_options;

namespace au = Acts::units;

namespace FW {

namespace Options {

  /// the particle gun options, the are prefixes with gp
  template <class AOPT>
  void
  addDD4hepOptions(AOPT& opt)
  {
    opt.add_options()(
        "dd4hep-input",
        po::value<read_strings>()->multitoken()->default_value(
            {"file:Detectors/DD4hepDetector/compact/FCChhBaseline/"
             "FCChh_DectEmptyMaster.xml",
             "file:Detectors/DD4hepDetector/compact/FCChhBaseline/"
             "FCChhTrackerTkLayout.xml"}),
        "The locations of the input DD4hep files, use 'file:foo.xml'. In case "
        "you want to read in multiple files, just seperate the strings by "
        "space.")("dd4hep-envelopeR",
                  po::value<double>()->default_value(0.),
                  "The envelop cover in R for DD4hep volumes.")(
        "dd4hep-envelopeR",
        po::value<double>()->default_value(1. * Acts::units::_mm),
        "The tolerance added to the geometrical extension in r of the "
        "layers contained to build the volume envelope around in mm.")(
        "dd4hep-envelopeZ",
        po::value<double>()->default_value(1. * Acts::units::_mm),
        "The tolerance added to the geometrical extension in z of the "
        "layers contained to build the volume envelope around in mm.")(
        "dd4hep-digitizationmodules",
        po::value<bool>()->default_value(false),
        "The envelop cover in z for DD4hep volumes.")(
        "dd4hep-loglevel",
        po::value<size_t>()->default_value(2),
        "The output log level of the geometry building. Please set the "
        "wished "
        "number (0 = VERBOSE, 1 = "
        "DEBUG, 2 = INFO, 3 = WARNING, 4 = ERROR, 5 = FATAL).");
  }

  /// read the particle gun options and return a Config file
  template <class AMAP>
  FW::DD4hep::GeometryService::Config
  readDD4hepConfig(const AMAP& vm)
  {
    Acts::Logging::Level logLevel
        = Acts::Logging::Level(vm["dd4hep-loglevel"].template as<size_t>());
    std::cout << "- the geometry building output log level is set to "
              << logLevel << std::endl;
    // DETECTOR configuration:
    // --------------------------------------------------------------------------------
    FW::DD4hep::GeometryService::Config gsConfig("GeometryService", logLevel);
    gsConfig.xmlFileNames = vm["dd4hep-input"].template as<read_strings>();
    gsConfig.bTypePhi     = Acts::equidistant;
    gsConfig.bTypeR       = Acts::arbitrary;
    gsConfig.bTypeZ       = Acts::equidistant;
    gsConfig.envelopeR    = vm["dd4hep-envelopeR"].template as<double>();
    gsConfig.envelopeZ    = vm["dd4hep-envelopeZ"].template as<double>();
    gsConfig.buildDigitizationModules
        = vm["dd4hep-digitizationmodules"].template as<bool>();
    return gsConfig;
  }
}
}

#endif  // ACTFW_OPTIONS_DD4HEPDETECTOROPTIONS_HPP

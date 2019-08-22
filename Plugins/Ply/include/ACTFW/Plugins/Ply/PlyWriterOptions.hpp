// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Geometry/TrackingGeometry.hpp>
#include <iostream>
#include "ACTFW/Plugins/Ply/PlyWriter.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace po = boost::program_options;

namespace FW {

namespace Options {

  /// Common ply writing options
  ///
  /// @tparam aopt_t Type of the options object (from BOOST)
  ///
  /// @param opt The options object, where string based options are attached
  template <typename aopt_t>
  void
  addPlyWriterOptions(aopt_t& opt)
  {
    opt.add_options()(
        "ply-output-file",
        po::value<std::string>()->default_value("DetectorMaterial.ply"),
        "The name of the output file.")(
        "ply-colour-range",
        po::value<std::vector<int>>()->multitoken()->default_value(
            std::vector<int>{0, 0, 0, 255, 255, 255}),
        "The range of colours to be used for rendering different materials "
        "(two RGB values).");
  }

  /// Read the program arguments and return a config object
  template <class AMAP>
  FW::Ply::PlyWriter::Config
  readPlyWriterConfig(const AMAP&                   vm,
                      const Acts::TrackingGeometry& tGeo,
                      const std::string&            lname = "PlyWriter",
                      Acts::Logging::Level loglevel       = Acts::Logging::INFO)
  {
    FW::Ply::PlyWriter::Config plyConfig(
        tGeo, "DetectorMaterial.ply", lname, loglevel);

    plyConfig.filename = vm["ply-output-file"].template as<std::string>();

    std::vector<int> cr
        = vm["ply-colour-range"].template as<std::vector<int>>();
    if (cr.size() != 6) cr = std::vector<int>{0, 0, 0, 255, 255, 255};

    plyConfig.minCol = std::array<int, 3>{{cr[0], cr[1], cr[2]}};
    plyConfig.maxCol = std::array<int, 3>{{cr[3], cr[4], cr[5]}};

    return plyConfig;
  }

}  // namespace Options
}  // namespace FW

// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>

#include <Acts/Material/IMaterialDecorator.hpp>
#include <Acts/Utilities/Logger.hpp>
#include "ACTFW/Geometry/MaterialWiper.hpp"
#include "ACTFW/Plugins/Json/JsonGeometryConverter.hpp"
#include "ACTFW/Plugins/Json/JsonMaterialDecorator.hpp"
#include "ACTFW/Plugins/Root/RootMaterialDecorator.hpp"

namespace FW {
namespace Geometry {

  /// @brief helper method to setup the geometry
  ///
  /// @tparam options_map_t Type of the options to be read
  /// @tparam geometry_setupt_t Type of the callable geometry setup
  ///
  /// @param vm the parsed options map
  /// @param geometrySetup the callable geometry setup
  ///
  /// @return a pair of TrackingGeometry and context decorators
  template <typename options_map_t, typename geometry_setup_t>
  auto
  build(const options_map_t& vm, geometry_setup_t& geometrySetup)
  {

    // Material decoration
    std::shared_ptr<const Acts::IMaterialDecorator> matDeco = nullptr;
    auto matType = vm["mat-input-type"].template as<std::string>();
    if (matType == "none") {
      matDeco = std::make_shared<const Acts::MaterialWiper>();
    } else if (matType == "file") {
      // Retrieve the filename
      auto fileName = vm["mat-input-file"].template as<std::string>();
      // json or root based decorator
      if (fileName.find(".json") != std::string::npos) {
        // Set up the converter first
        FW::Json::JsonGeometryConverter::Config jsonGeoConvConfig;
        // Set up the json-based decorator
        matDeco = std::make_shared<const FW::Json::JsonMaterialDecorator>(
            jsonGeoConvConfig, fileName);
      } else if (fileName.find(".root") != std::string::npos) {
        // Set up the root-based decorator
        FW::Root::RootMaterialDecorator::Config rootMatDecConfig;
        rootMatDecConfig.fileName = fileName;
        matDeco = std::make_shared<const FW::Root::RootMaterialDecorator>(
            rootMatDecConfig);
      }
    }

    /// Return the geometry and context decorators
    return geometrySetup(vm, matDeco);
  }

}  // namespace Geometry
}  // namespace FW

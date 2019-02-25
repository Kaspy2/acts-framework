// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/GenericDetector/BuildGenericDetector.hpp"
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "detail/VoidMaterialReader.hpp"

/// @brief adding some specific options for this geometry type
struct GenericOptions
{
  /// @brief operator to be called to add options for the generic detector
  ///
  // @tparam options_t Type of the options object
  ///@param opt Options object
  template <typename options_t>
  void
  operator()(options_t& opt)
  {
  }
};

/// @brief geometry getter, the operator() will be called int he example base
struct GenericGeometry
{
  /// @brief operator called to construct the tracking geometry
  ///
  /// @tparam variable_map_t Type of the variable map template for parameters
  /// @tparam material_reader_t the source for the surface material map
  ///
  /// @param vm the parameter map object
  ///
  /// @return a closed TrackingGeometry object
  template <typename variable_map_t,
            typename material_reader_t = VoidMaterialReader>
  std::shared_ptr<const Acts::TrackingGeometry>
  operator()(variable_map_t&          vm,
             const material_reader_t& mreader = VoidMaterialReader())
  {
    // --------------------------------------------------------------------------------
    // Set geometry building logging level
    Acts::Logging::Level surfaceLogLevel = Acts::Logging::Level(
        vm["geo-surface-loglevel"].template as<size_t>());
    Acts::Logging::Level layerLogLevel
        = Acts::Logging::Level(vm["geo-layer-loglevel"].template as<size_t>());
    Acts::Logging::Level volumeLogLevel
        = Acts::Logging::Level(vm["geo-volume-loglevel"].template as<size_t>());

    // Check if material should be built
    size_t materialMode = vm["geo-material-mode"].template as<size_t>();

    /// return the generic detector
    return FW::Generic::buildGenericDetector(int(materialMode),
                                             mreader(),
                                             surfaceLogLevel,
                                             layerLogLevel,
                                             volumeLogLevel,
                                             3);
  }
};
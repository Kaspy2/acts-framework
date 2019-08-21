// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <fstream>
#include <map>
#include <mutex>
#include "ACTFW/Plugins/Json/JsonGeometryConverter.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"

// Convenience shorthand
using json = nlohmann::json;

namespace Acts {

using SurfaceMaterialMap
    = std::map<GeometryID, std::shared_ptr<const ISurfaceMaterial>>;

using VolumeMaterialMap
    = std::map<GeometryID, std::shared_ptr<const IVolumeMaterial>>;
}  // namespace Acts

namespace FW {

namespace Json {

  /// @brief Material decorator from Json format
  ///
  /// This reads in material maps for surfaces and volumes
  /// from a json file
  class JsonMaterialDecorator : public Acts::IMaterialDecorator
  {
  public:
    JsonMaterialDecorator(const JsonGeometryConverter::Config& rConfig,
                          const std::string&                   jFileName,
                          bool clearSurfaceMaterial = true,
                          bool clearVolumeMaterial  = true)
      : m_readerConfig(rConfig)
      , m_clearSurfaceMaterial(clearSurfaceMaterial)
      , m_clearVolumeMaterial(clearVolumeMaterial)
    {
      // the material reader
      FW::Json::JsonGeometryConverter::Config jmConverterCfg(
          "JsonGeometryConverter", Acts::Logging::VERBOSE);
      FW::Json::JsonGeometryConverter jmConverter(jmConverterCfg);

      std::ifstream ifj(jFileName.c_str());
      json          jin;
      ifj >> jin;

      auto maps            = jmConverter.jsonToMaterialMaps(jin);
      m_surfaceMaterialMap = maps.first;
      m_volumeMaterialMap  = maps.second;
    }

    /// Decorate a surface
    ///
    /// @param surface the non-cost surface that is decorated
    void
    decorate(Acts::Surface& surface) const final
    {
      // Clear the material if registered to do so
      if (m_clearSurfaceMaterial) { surface.assignSurfaceMaterial(nullptr); }

      // Try to find the surface in the map
      auto sMaterial = m_surfaceMaterialMap.find(surface.geoID());
      if (sMaterial != m_surfaceMaterialMap.end()) {
        surface.assignSurfaceMaterial(sMaterial->second);
      }
    }

    /// Decorate a TrackingVolume
    ///
    /// @param volume the non-cost volume that is decorated
    void
    decorate(Acts::TrackingVolume& volume) const final
    {
      // Clear the material if registered to do so
      if (m_clearVolumeMaterial) { volume.assignVolumeMaterial(nullptr); }
    }

  private:
    JsonGeometryConverter::Config m_readerConfig;
    Acts::SurfaceMaterialMap      m_surfaceMaterialMap;
    Acts::VolumeMaterialMap       m_volumeMaterialMap;

    bool m_clearSurfaceMaterial{true};
    bool m_clearVolumeMaterial{true};
  };

}  // namespace Json
}  // namespace FW

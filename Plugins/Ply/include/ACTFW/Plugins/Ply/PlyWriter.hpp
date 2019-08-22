// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <fstream>
#include <vector>
#include "ACTFW/Framework/AlgorithmContext.hpp"
#include "ACTFW/Framework/IWriter.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/IVisualization.hpp"
#include "Acts/Utilities/PlyHelper.hpp"

namespace Acts {
class TrackingVolume;
class TrackingGeometry;
}  // namespace Acts

namespace FW {

namespace Ply {

  class PlyWriter : public FW::IWriter
  {

  public:
    Acts::PlyHelper<double> m_plyh;
    std::string             m_filename;
    float                   m_min;
    float                   m_max;

    // @class Config
    //
    // The nested config class
    class Config
    {
    public:
      /// the default logger
      std::shared_ptr<const Acts::Logger> logger;
      /// the name of the writer
      std::string name = "";

      /// TrackingGeometry object reference
      const Acts::TrackingGeometry* tGeometry;
      /// output directory for the .ply file
      std::string filename = "";

      /// minimum and maximum colours (colour range for material rendering)
      Acts::IVisualization::color_type minCol = {0, 0, 0};
      Acts::IVisualization::color_type maxCol = {255, 255, 255};

      Config(const Acts::TrackingGeometry& tGeo,
             const std::string             fn    = "DetectorMaterial.ply",
             const std::string&            lname = "PlyWriter",
             Acts::Logging::Level          lvl   = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, lvl)), name(lname), filename(fn)
      {
        tGeometry = &tGeo;
      };
    };

    /// Constructor
    /// @param cfg is the configuration class
    PlyWriter(const Config& cfg) : m_cfg(cfg), m_filename(cfg.filename)
    {
      m_plyh = Acts::PlyHelper<double>();
    };

    std::string
    name() const override;

    FW::ProcessCode
    write(Acts::GeometryContext         context,
          const Acts::TrackingGeometry& tGeometry);

    void
    write(Acts::GeometryContext context, const Acts::TrackingVolume& tVolume);

    FW::ProcessCode
    write(Acts::GeometryContext context, const Acts::Surface& surface);

    FW::ProcessCode
    finalize(std::ostream& os);

    Acts::IVisualization::color_type
    getColor(const Acts::ISurfaceMaterial* sMaterial);

    Acts::IVisualization::color_type
    getColor(float val);

    void
    writeCylinder(unsigned int                     nSegments,
                  const Acts::Transform3D&         transform,
                  double                           r,
                  double                           halfZ,
                  Acts::IVisualization::color_type color);

    void
    writeDisk(unsigned int                     nSegments,
              const Acts::Transform3D&         transform,
              double                           rMin,
              double                           rMax,
              Acts::IVisualization::color_type color);

    void
    writeCylinderBinned(const Acts::Transform3D&           transform,
                        double                             r,
                        double                             halfZ,
                        const Acts::BinnedSurfaceMaterial& bsm);

    void
    writeDiskBinned(const Acts::Transform3D&           transform,
                    double                             rMin,
                    double                             rMax,
                    const Acts::BinnedSurfaceMaterial& bsm);

    void
    writeBinnedSurface(const Acts::Surface& surface);

    FW::ProcessCode
    write(const AlgorithmContext& context) final override;

    FW::ProcessCode
    endRun() final override;

  private:
    Config m_cfg;  ///< the config class

    /// Private access to the logging instance
    const Acts::Logger&
    logger() const
    {
      return *m_cfg.logger;
    }

    FW::ProcessCode
    getMinMax(Acts::GeometryContext         context,
              const Acts::TrackingGeometry& tGeometry,
              float&                        upper,
              float&                        lower);

    void
    getMinMax(Acts::GeometryContext       context,
              const Acts::TrackingVolume& tVolume,
              float&                      upper,
              float&                      lower);

    FW::ProcessCode
    getMinMax(Acts::GeometryContext context,
              const Acts::Surface&  surface,
              float&                upper,
              float&                lower);
  };

}  // namespace Ply
}  // namespace FW

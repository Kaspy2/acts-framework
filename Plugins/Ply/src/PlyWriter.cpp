// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Ply/PlyWriter.hpp"

#include <cmath>
#include <ios>
#include <iostream>
#include <limits>
#include <stdexcept>

#include <Acts/Geometry/Layer.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <Acts/Geometry/GeometryID.hpp>
#include <Acts/Surfaces/CylinderBounds.hpp>
#include <Acts/Surfaces/PlanarBounds.hpp>
#include <Acts/Surfaces/RadialBounds.hpp>
#include <Acts/Surfaces/SurfaceBounds.hpp>

FW::ProcessCode
FW::Ply::PlyWriter::finalize(std::ostream& os)
{
  m_plyh.write(os);
  return FW::ProcessCode::SUCCESS;
};

// process the tracking geometry
FW::ProcessCode
FW::Ply::PlyWriter::write(Acts::GeometryContext         context,
                          const Acts::TrackingGeometry& tGeometry)
{
  ACTS_DEBUG(">>Ply: Writer for TrackingGeometry object called.");
  // get the world volume
  auto world = tGeometry.highestTrackingVolume();
  if (world) write(context, *world);
  // return the success code
  return FW::ProcessCode::SUCCESS;
}

/// process this volume
void
FW::Ply::PlyWriter::write(Acts::GeometryContext       context,
                          const Acts::TrackingVolume& tVolume)
{
  ACTS_DEBUG(">>Ply: Writer for TrackingVolume object called.");
  // get the confined layers and process them
  if (tVolume.confinedLayers()) {
    ACTS_VERBOSE(">>Ply: Layers are present, process them.");
    // loop over the layers
    for (auto layer : tVolume.confinedLayers()->arrayObjects()) {
      // we jump navigation layers
      if (layer->layerType() == Acts::navigation) continue;

      // try to write the material surface
      if (layer->surfaceRepresentation().surfaceMaterial()) {
        write(context, layer->surfaceRepresentation());
      }

      // the approaching surfaces and check if they have material
      if (layer->approachDescriptor()) {
        // loop over the contained Surfaces
        for (auto& cSurface : layer->approachDescriptor()->containedSurfaces())
          if (cSurface->surfaceMaterial()) {
            write(context, *cSurface);
          }
      }
      // check for sensitive surfaces
      if (layer->surfaceArray()) {
        ACTS_VERBOSE(">>Ply: There are "
                     << layer->surfaceArray()->surfaces().size()
                     << " surfaces.");

        // loop over the surfaces
        for (auto& surface : layer->surfaceArray()->surfaces()) {
          if (surface && (write(context, *surface)) == FW::ProcessCode::ABORT)
            return;
        }
      }
    }
  }
  // Recursive self call
  // get the confined volumes and step down the hierarchy
  if (tVolume.confinedVolumes()) {
    // loop over the volumes and write what they have
    for (auto volume : tVolume.confinedVolumes()->arrayObjects()) {
      write(context, *volume.get());
    }
  }
}

/// process this surface
FW::ProcessCode
FW::Ply::PlyWriter::write(Acts::GeometryContext context,
                          const Acts::Surface&  surface)
{
  ACTS_DEBUG(">>Ply: Writer for Surface object called.");

  bool outputSensitive    = true;
  bool outputLayerSurface = true;
  bool outputThickness    = false;

  bool foundBinnedSurface = false;

  // let's get the bounds & the transform
  const Acts::SurfaceBounds& surfaceBounds = surface.bounds();
  auto                       sTransform    = surface.transform(context);

  auto sMaterial = surface.surfaceMaterial();

  if (sMaterial != nullptr) {
    /// Found surface material

    auto bsMaterial
        = dynamic_cast<const Acts::BinnedSurfaceMaterial*>(sMaterial);
    if (bsMaterial != nullptr) {
      // Found a binned surface
      foundBinnedSurface = true;
    }

    // dynamic_cast to PlanarBounds
    const Acts::PlanarBounds* planarBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&surfaceBounds);
    // only continue if the cast worked
    if (planarBounds && outputSensitive) {
      ACTS_DEBUG(">>Ply: Writing out a PlaneSurface");
      // get the vertices
      auto planarVertices = planarBounds->vertices();
      // loop over the vertices
      std::vector<Acts::Vector3D> vertices;
      vertices.reserve(planarVertices.size());
      for (auto pv : planarVertices) {
        // get the point in 3D
        Acts::Vector3D v3D(sTransform * Acts::Vector3D(pv.x(), pv.y(), 0.));
        vertices.push_back(v3D);
      }
      // output to file
      m_plyh.face(vertices, getColor(sMaterial));
    }

    // check if you have layer and check what you have
    // dynamic cast to CylinderBounds work the same
    const Acts::CylinderBounds* cylinderBounds
        = dynamic_cast<const Acts::CylinderBounds*>(&surfaceBounds);
    if (cylinderBounds && outputLayerSurface) {
      if (foundBinnedSurface) {
        ACTS_DEBUG(">>Ply: Writing out a binned CylinderSurface");
        writeCylinderBinned(sTransform,
                            cylinderBounds->r(),
                            cylinderBounds->halflengthZ(),
                            *bsMaterial);
        // writeCylinder(100, sTransform, cylinderBounds->r(),
        // cylinderBounds->halflengthZ(), getColor(sMaterial));
      } else {
        ACTS_DEBUG(">>Ply: Writing out a CylinderSurface");
        writeCylinder(100,
                      sTransform,
                      cylinderBounds->r(),
                      cylinderBounds->halflengthZ(),
                      getColor(sMaterial));
      }
    }

    // dynamic cast to RadialBounds or disc bounds work the same
    const Acts::RadialBounds* radialBounds
        = dynamic_cast<const Acts::RadialBounds*>(&surfaceBounds);
    if (radialBounds && outputLayerSurface) {
      if (foundBinnedSurface) {
        ACTS_DEBUG(">>Ply: Writing out a binned DiskSurface");
        writeDiskBinned(sTransform,
                        radialBounds->rMin(),
                        radialBounds->rMax(),
                        *bsMaterial);
      } else {
        ACTS_DEBUG(">>Ply: Writing out a DiskSurface");
        writeDisk(100,
                  sTransform,
                  radialBounds->rMin(),
                  radialBounds->rMax(),
                  getColor(sMaterial));
      }
    }

    // if (foundBinnedSurface) {
    //   writeBinnedSurface(surface);
    // }

  } else {
    /// No surface material
    /// Could still draw but with a fixed colour
  }

  // return success
  return FW::ProcessCode::SUCCESS;
}

/// translate surface material into colour by type (homogeneous, proto, binned)
Acts::IVisualization::color_type
FW::Ply::PlyWriter::getColor(const Acts::ISurfaceMaterial* sMaterial)
{
  auto psMaterial = dynamic_cast<const Acts::ProtoSurfaceMaterial*>(sMaterial);
  if (psMaterial != nullptr) {
    // type is proto
    return {0, 255, 0};
  }

  auto hsMaterial
      = dynamic_cast<const Acts::HomogeneousSurfaceMaterial*>(sMaterial);
  if (hsMaterial != nullptr) {
    // type is homogeneous
    return {255, 0, 0};
  }

  auto bsMaterial = dynamic_cast<const Acts::BinnedSurfaceMaterial*>(sMaterial);
  if (bsMaterial != nullptr) {
    // type is binned
    return {0, 0, 255};
  }

  return {0, 0, 0};
}

/// calculate colour from surface material's radiation length X0
Acts::IVisualization::color_type
FW::Ply::PlyWriter::getColor(float val)
{

  if (val > m_max || val < m_min)
    return {0, 0, 0};
  else {
    // a simple linear interpolation between min and max colour
    float normedVal = std::log((val - m_min)) / std::log((m_max - m_min));
    return {
        static_cast<int>(normedVal * m_cfg.maxCol[0]
                         + (1 - normedVal) * m_cfg.minCol[0]),
        static_cast<int>(normedVal * m_cfg.maxCol[1]
                         + (1 - normedVal) * m_cfg.minCol[1]),
        static_cast<int>(normedVal * m_cfg.maxCol[2]
                         + (1 - normedVal) * m_cfg.minCol[2]),
    };
  }
  // else {
  //   float r = 0. , g = 0. , b = 0. ;

  //   float redLim = 0.4;
  //   float bluLim = 0.6;
  //   float lgrLim = 0.3;
  //   float ugrLim = 0.7;

  //   // normalize val to a value between 0 and 1
  //   float normedVal = (val - m_min) / (m_max - m_min);

  //   if (normedVal < redLim) {
  //     r = m_cfg.minCol[0] +
  //     (m_cfg.maxCol[0]-m_cfg.minCol[0])*(normedVal-0.)/(redLim-0.);
  //   }
  //   if (normedVal > bluLim) {
  //     b = m_cfg.minCol[2] +
  //     (m_cfg.maxCol[2]-m_cfg.minCol[2])*(normedVal-bluLim)/(1.-bluLim);
  //   }
  //   if ((normedVal > lgrLim) && (normedVal < ugrLim)) {
  //     g = m_cfg.minCol[1] +
  //     (m_cfg.maxCol[1]-m_cfg.minCol[1])*(normedVal-lgrLim)/(ugrLim-lgrLim);
  //   }

  //   return {static_cast<int>(r), static_cast<int>(g), static_cast<int>(b)};

  // }
}

void
FW::Ply::PlyWriter::writeCylinder(unsigned int                     nSegments,
                                  const Acts::Transform3D&         transform,
                                  double                           r,
                                  double                           halfZ,
                                  Acts::IVisualization::color_type color)
{
  double phistep = 2 * M_PI / nSegments;

  Acts::Vector3D point1(transform
                        * Acts::Vector3D(r * cos(0), r * sin(0), halfZ));
  Acts::Vector3D point2(transform
                        * Acts::Vector3D(r * cos(0), r * sin(0), -halfZ));

  // loop over phi steps
  for (size_t iphi = 1; iphi <= nSegments; ++iphi) {
    double phi = iphi * phistep;

    Acts::Vector3D point3(transform
                          * Acts::Vector3D(r * cos(phi), r * sin(phi), halfZ));
    Acts::Vector3D point4(transform
                          * Acts::Vector3D(r * cos(phi), r * sin(phi), -halfZ));

    std::vector<Acts::Vector3D> vertices = {point1, point2, point4, point3};

    // output to file
    m_plyh.face(vertices, color);

    point1 = point3;
    point2 = point4;
  }
}

void
FW::Ply::PlyWriter::writeDisk(unsigned int                     nSegments,
                              const Acts::Transform3D&         transform,
                              double                           rMin,
                              double                           rMax,
                              Acts::IVisualization::color_type color)
{
  double phistep = 2 * M_PI / nSegments;

  Acts::Vector3D point1(transform
                        * Acts::Vector3D(rMin * cos(0), rMin * sin(0), 0));
  Acts::Vector3D point2(transform
                        * Acts::Vector3D(rMax * cos(0), rMax * sin(0), 0));

  // loop over phi steps
  for (size_t iphi = 1; iphi <= nSegments; ++iphi) {
    double phi = iphi * phistep;

    Acts::Vector3D point3(
        transform * Acts::Vector3D(rMin * cos(phi), rMin * sin(phi), 0));
    Acts::Vector3D point4(
        transform * Acts::Vector3D(rMax * cos(phi), rMax * sin(phi), 0));

    std::vector<Acts::Vector3D> vertices = {point1, point2, point4, point3};

    // output to file
    m_plyh.face(vertices, color);

    point1 = point3;
    point2 = point4;
  }
}

void
FW::Ply::PlyWriter::writeCylinderBinned(const Acts::Transform3D& transform,
                                        double                   r,
                                        double                   halfZ,
                                        const Acts::BinnedSurfaceMaterial& bsm)
{

  auto binUtil = bsm.binUtility();

  int phiSegments = binUtil.bins(0);
  int zSegments   = binUtil.bins(1);

  // std::cout<< ">>Ply:BinPhi:"<<
  // Acts::binningValueNames[binUtil.binningValue(0)] <<std::endl;
  // std::cout<< ">>Ply:BinZ  :"<<
  // Acts::binningValueNames[binUtil.binningValue(1)] <<std::endl;

  double phistep = 2 * M_PI / phiSegments;
  double zstep   = 2 * halfZ / zSegments;

  std::vector<Acts::Vector3D> phiPoints;
  phiPoints.reserve(zSegments + 1);

  for (int iz = 0; iz <= zSegments; iz++) {
    double         z = -halfZ + (iz * zstep);
    Acts::Vector3D point(transform * Acts::Vector3D(r * cos(0), r * sin(0), z));
    phiPoints.emplace_back(point);
  }

  // loop over phi steps
  for (int iphi = 1; iphi <= phiSegments; iphi++) {
    double phi = iphi * phistep;

    Acts::Vector3D point3(transform
                          * Acts::Vector3D(r * cos(phi), r * sin(phi), -halfZ));

    for (int iz = 1; iz <= zSegments; iz++) {
      double z = -halfZ + (iz * zstep);

      Acts::Vector3D point4(transform
                            * Acts::Vector3D(r * cos(phi), r * sin(phi), z));

      std::vector<Acts::Vector3D> vertices
          = {phiPoints[iz - 1], phiPoints[iz], point4, point3};

      auto col = bsm.materialProperties(iphi - 1, iz - 1).averageX0();
      // auto col = bsm.materialProperties(iphi-1, iz-1).averageA();

      Acts::IVisualization::color_type color = getColor(col);
      // Acts::IVisualization::color_type color = {0, 255, 0};
      // if (col == 0){
      //   color = {0, 0, 0};
      // }

      // output to file
      m_plyh.face(vertices, color);

      phiPoints[iz - 1] = point3;
      point3            = point4;
    }

    phiPoints[zSegments] = point3;
  }
}

void
FW::Ply::PlyWriter::writeDiskBinned(const Acts::Transform3D& transform,
                                    double                   rMin,
                                    double                   rMax,
                                    const Acts::BinnedSurfaceMaterial& bsm)
{

  auto binUtil = bsm.binUtility();

  int rSegments   = binUtil.bins(0);
  int phiSegments = binUtil.bins(1);

  // std::cout<< ">>Ply:BinR  :"<<
  // Acts::binningValueNames[binUtil.binningValue(0)] <<std::endl;
  // std::cout<< ">>Ply:BinPhi:"<<
  // Acts::binningValueNames[binUtil.binningValue(1)] <<std::endl;

  double phistep = 2 * M_PI / phiSegments;
  double rstep   = (rMax - rMin) / rSegments;

  std::vector<Acts::Vector3D> phiPoints;
  phiPoints.reserve(rSegments + 1);

  for (int ir = 0; ir <= rSegments; ir++) {
    double         r = rMin + (ir * rstep);
    Acts::Vector3D point(transform * Acts::Vector3D(r * cos(0), r * sin(0), 0));
    phiPoints.emplace_back(point);
  }

  // loop over phi steps
  for (int iphi = 1; iphi <= phiSegments; iphi++) {
    double phi = iphi * phistep;

    Acts::Vector3D point3(
        transform * Acts::Vector3D(rMin * cos(phi), rMin * sin(phi), 0));

    for (int ir = 1; ir <= rSegments; ir++) {
      double r = rMin + (ir * rstep);

      Acts::Vector3D point4(transform
                            * Acts::Vector3D(r * cos(phi), r * sin(phi), 0));

      std::vector<Acts::Vector3D> vertices
          = {phiPoints[ir - 1], phiPoints[ir], point4, point3};

      auto col = bsm.materialProperties(ir - 1, iphi - 1).averageX0();
      // auto col = bsm.materialProperties(ir-1, iphi-1).averageA();

      Acts::IVisualization::color_type color = getColor(col);
      // std::array<int,3> color = {0, 255, 0};
      // if (col == 0){
      //   color = {0, 0, 0};
      // }

      // output to file
      m_plyh.face(vertices, color);

      phiPoints[ir - 1] = point3;
      point3            = point4;
    }

    phiPoints[rSegments] = point3;
  }
}

/// An attempt at creating a general surface writer
/// (failed since equidistant binning does not define boundaries in BinUtility).
void
FW::Ply::PlyWriter::writeBinnedSurface(const Acts::Surface& surface)
{

  const Acts::BinnedSurfaceMaterial* bSurfaceMaterial
      = dynamic_cast<const Acts::BinnedSurfaceMaterial*>(
          surface.surfaceMaterial());
  const Acts::SurfaceBounds& bSurfaceBounds = surface.bounds();
  const Acts::BinUtility&    binUtil        = bSurfaceMaterial->binUtility();

  int xSegments = binUtil.bins(0);
  int ySegments = binUtil.bins(1);
  int zSegments = binUtil.bins(2);  // optional

  // Acts::Vector3D lpos(0,0,0);
  // Acts::Vector3D gpos;
  // surface.localToGlobal(0, lpos, 0, gpos);

  auto binningData = binUtil.binningData();

  auto binningDataX = binningData[0];
  auto binningDataY = binningData[1];

  auto binBoundsX = binningDataX.boundaries();
  auto binBoundsY = binningDataY.boundaries();

  std::cout << binUtil << std::endl;  // binUtility must have been called with
                                      // the basic equidistant constructor

  for (auto v : bSurfaceBounds.valueStore()) {
    std::cout << v << " ";
  }

  std::cout << std::endl;

  // no info on boundaries or step size or min/max

  // if (binningDataX.type == Acts::BinningType::equidistant){
  //   std::cout<< "PlyGREP: X: Equidistant , ";
  // }
  // else{
  //   std::cout<< "PlyGREP: X: Arbitrary , ";
  // }
  // if (binningDataX.option == Acts::BinningOption::closed) {
  //   std::cout<< "Closed" << std::endl;
  // }
  // else{
  //   std::cout<< "Open" << std::endl;
  // }

  // if (binningDataY.type == Acts::BinningType::equidistant){
  //   std::cout<< "PlyGREP: Y: Equidistant , ";
  // }
  // else{
  //   std::cout<< "PlyGREP: Y: Arbitrary , ";
  // }
  // if (binningDataY.option == Acts::BinningOption::closed) {
  //   std::cout<< "Closed" << std::endl;
  // }
  // else{
  //   std::cout<< "Open" << std::endl;
  // }

  // std::cout<<"PlyGREP: "<<Acts::binningValueNames[binUtil.binningValue(0)]<<"
  // : ";

  // for (auto bx: binBoundsX) {
  //   std::cout<<bx<<" ";
  // }
  // std::cout<<std::endl;

  // std::cout<<"PlyGREP: "<<Acts::binningValueNames[binUtil.binningValue(1)]<<"
  // : ";

  // for (auto by: binBoundsY) {
  //   std::cout<<by<<" ";
  // }
  // std::cout<<std::endl;
}

std::string
FW::Ply::PlyWriter::name() const
{
  return m_cfg.name;
};

/// Per-event writer (do nothing).
FW::ProcessCode
FW::Ply::PlyWriter::write(const AlgorithmContext& context)
{
  return FW::ProcessCode::SUCCESS;
}

/// Main runner (get material X0 bounds, recursively explore the geometry,
/// write out cylinders/bins with colour dependent on material X0 per bin,
/// output to file.)
FW::ProcessCode
FW::Ply::PlyWriter::endRun()
{
  float upper = -std::numeric_limits<float>::infinity();
  float lower = std::numeric_limits<float>::infinity();

  getMinMax(0, *(m_cfg.tGeometry), upper, lower);

  m_min = lower;
  m_max = upper;

  write(0, *(m_cfg.tGeometry));

  ACTS_INFO("X0 bounds : (" << lower << "," << upper << ")");

  ACTS_INFO("Saved in " << m_cfg.filename);

  auto outputStream = std::shared_ptr<std::ofstream>(new std::ofstream);
  outputStream->open(m_cfg.filename, std::ofstream::out | std::ofstream::trunc);
  finalize(*outputStream);
  outputStream->close();

  return FW::ProcessCode::SUCCESS;
};

// process the tracking geometry
FW::ProcessCode
FW::Ply::PlyWriter::getMinMax(Acts::GeometryContext         context,
                              const Acts::TrackingGeometry& tGeometry,
                              float&                        upper,
                              float&                        lower)
{
  // get the world volume
  auto world = tGeometry.highestTrackingVolume();
  if (world) getMinMax(context, *world, upper, lower);
  // return the success code
  return FW::ProcessCode::SUCCESS;
}

/// process this volume
void
FW::Ply::PlyWriter::getMinMax(Acts::GeometryContext       context,
                              const Acts::TrackingVolume& tVolume,
                              float&                      upper,
                              float&                      lower)
{
  // get the confined layers and process them
  if (tVolume.confinedLayers()) {
    // loop over the layers
    for (auto layer : tVolume.confinedLayers()->arrayObjects()) {
      // we jump navigation layers
      if (layer->layerType() == Acts::navigation) continue;

      // explore the material surface
      if (layer->surfaceRepresentation().surfaceMaterial()) {
        getMinMax(context, layer->surfaceRepresentation(), upper, lower);
      }

      // the approaching surfaces and check if they have material
      if (layer->approachDescriptor()) {
        // loop over the contained Surfaces
        for (auto& cSurface : layer->approachDescriptor()->containedSurfaces())
          if (cSurface->surfaceMaterial()) {
            getMinMax(context, *cSurface, upper, lower);
          }
      }
      // check for sensitive surfaces
      if (layer->surfaceArray()) {
        // loop over the surfaces
        for (auto& surface : layer->surfaceArray()->surfaces()) {
          if (surface
              && (getMinMax(context, *surface, upper, lower))
                  == FW::ProcessCode::ABORT)
            return;
        }
      }
    }
  }
  // Recursive self call
  // get the confined volumes and step down the hierarchy
  if (tVolume.confinedVolumes()) {
    // loop over the volumes
    for (auto volume : tVolume.confinedVolumes()->arrayObjects()) {
      getMinMax(context, *volume.get(), upper, lower);
    }
  }
}

/// process this surface
FW::ProcessCode
FW::Ply::PlyWriter::getMinMax(Acts::GeometryContext context,
                              const Acts::Surface&  surface,
                              float&                upper,
                              float&                lower)
{
  const Acts::SurfaceBounds& surfaceBounds = surface.bounds();
  auto                       sTransform    = surface.transform(context);

  auto sMaterial = surface.surfaceMaterial();

  if (sMaterial != nullptr) {
    auto bsMaterial
        = dynamic_cast<const Acts::BinnedSurfaceMaterial*>(sMaterial);
    if (bsMaterial != nullptr) {
      auto butil = bsMaterial->binUtility();

      for (size_t i = 0; i < butil.bins(0); i++) {
        for (size_t j = 0; j < butil.bins(1); j++) {

          float curval = bsMaterial->materialProperties(i, j).averageX0();
          if ((curval > upper)
              && (curval != std::numeric_limits<float>::infinity())) {
            upper = curval;
          }
          if ((curval < lower)
              && (curval != -std::numeric_limits<float>::infinity())) {
            lower = curval;
          }
        }
      }
    }
  }

  // return success
  return FW::ProcessCode::SUCCESS;
}
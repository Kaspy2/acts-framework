// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/program_options.hpp>
#include "ACTFW/Plugins/DD4hep/DD4hepDetectorOptions.hpp"
#include "ACTFW/Plugins/DD4hep/GeometryService.hpp"
#include "ACTFW/Plugins/Obj/ObjSurfaceWriter.hpp"
#include "ACTFW/Plugins/Obj/ObjTrackingGeometryWriter.hpp"
#include "Acts/Detector/TrackingGeometry.hpp"

namespace po = boost::program_options;

int
main(int argc, char* argv[])
{

  // Declare the supported program options.
  po::options_description desc("Allowed options");
  // add help message
  desc.add_options()("help", "Produce help message");
  // add the detector options
  FW::Options::addDD4hepOptions<po::options_description>(desc);
  po::variables_map vm;
  // Get all options from contain line and store it into the map
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  // print help if requested
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }
  // read the detector config & dd4hep detector
  auto dd4HepDetectorConfig
      = FW::Options::readDD4hepConfig<po::variables_map>(vm);
  auto geometrySvc
      = std::make_shared<FW::DD4hep::GeometryService>(dd4HepDetectorConfig);
  std::shared_ptr<const Acts::TrackingGeometry> dd4tGeometry
      = geometrySvc->trackingGeometry();

  // the detectors
  std::vector<std::string> subDetectors = {"beampipe",
                                           "FCChhInner0",
                                           "FCChhInner",
                                           "FCChhOuter",
                                           "FCChhForwardHelper",
                                           "FCChhForward"};
  // the writers
  std::vector<std::shared_ptr<FWObj::ObjSurfaceWriter>> subWriters;
  std::vector<std::shared_ptr<std::ofstream>>           subStreams;
  // loop and create
  for (auto sdet : subDetectors) {
    // sub detector stream
    auto        sdStream = std::shared_ptr<std::ofstream>(new std::ofstream);
    std::string sdOutputName = sdet + std::string(".obj");
    sdStream->open(sdOutputName);
    // object surface writers
    FWObj::ObjSurfaceWriter::Config sdObjWriterConfig(sdet,
                                                      Acts::Logging::INFO);
    sdObjWriterConfig.filePrefix         = "mtllib materials.mtl'\n'";
    sdObjWriterConfig.outputPhiSegemnts  = 72;
    sdObjWriterConfig.outputPrecision    = 6;
    sdObjWriterConfig.outputScalor       = 1.;
    sdObjWriterConfig.outputThickness    = 1.;
    sdObjWriterConfig.outputSensitive    = true;
    sdObjWriterConfig.outputLayerSurface = true;
    sdObjWriterConfig.outputStream       = sdStream;
    auto sdObjWriter
        = std::make_shared<FWObj::ObjSurfaceWriter>(sdObjWriterConfig);
    // push back
    subWriters.push_back(sdObjWriter);
    subStreams.push_back(sdStream);
  }
  // configure the tracking geometry writer
  FWObj::ObjTrackingGeometryWriter::Config tgObjWriterConfig(
      "ObjTrackingGeometryWriter", Acts::Logging::VERBOSE);
  tgObjWriterConfig.surfaceWriters       = subWriters;
  tgObjWriterConfig.filePrefix           = "mtllib materials.mtl'\n'";
  tgObjWriterConfig.sensitiveGroupPrefix = "usemtl silicon'\n'";
  tgObjWriterConfig.layerPrefix          = "usemtl support'\n'";
  // the tracking geometry writers
  auto tgObjWriter
      = std::make_shared<FWObj::ObjTrackingGeometryWriter>(tgObjWriterConfig);

  // write the tracking geometry object
  tgObjWriter->write(*(dd4tGeometry.get()));

  // --------------------------------------------------------------------------------
  // close the output streams
  for (auto sStreams : subStreams) {
    sStreams->close();
  }
}

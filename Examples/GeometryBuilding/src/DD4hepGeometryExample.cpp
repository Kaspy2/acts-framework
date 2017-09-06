///////////////////////////////////////////////////////////////////
// GenatinoRecording.cpp
///////////////////////////////////////////////////////////////////

#include <boost/program_options.hpp>
#include "ACTFW/Plugins/DD4hep/GeometryService.hpp"
#include "ACTFW/Plugins/Obj/ObjSurfaceWriter.hpp"
#include "ACTFW/Plugins/Obj/ObjTrackingGeometryWriter.hpp"
#include "ACTS/Detector/TrackingGeometry.hpp"

namespace po = boost::program_options;

int
main(int argc, char* argv[])
{
  
  // Declare the supported program options.
  po::options_description desc("Allowed options");
  desc.add_options()("help", "Produce help message")(
      "input",
      po::value<std::string>()->default_value(
        "file:Detectors/DD4hepDetector/compact/FCChhTrackerTkLayout.xml"),
      "The location of the input DD4hep file, use 'file:foo.xml'");  

  po::variables_map vm;
  // Get all options from contain line and store it into the map
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  // print help if needed
  // output messages
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }
  
  // DETECTOR:
  // --------------------------------------------------------------------------------
  // DD4Hep detector definition
  //
  // set up the geometry service
  FW::DD4hep::GeometryService::Config gsConfig("GeometryService",
                                             Acts::Logging::INFO);

  gsConfig.xmlFileName = vm["input"].as<std::string>();
  gsConfig.bTypePhi  = Acts::equidistant;
  gsConfig.bTypeR    = Acts::equidistant;
  gsConfig.bTypeZ    = Acts::equidistant;
  gsConfig.envelopeR = 0.;
  gsConfig.envelopeZ = 0.;

  auto geometrySvc = std::make_shared<FW::DD4hep::GeometryService>(gsConfig);
  std::shared_ptr<const Acts::TrackingGeometry> dd4Geometry
      = geometrySvc->trackingGeometry();

  // the detectors
  std::vector<std::string> subDetectors
      = {"beampipe", "FCChhInner0", "FCChhInner", "FCChhOuter", "FCChhForwardHelper", "FCChhForward"};
  // the writers
  std::vector<std::shared_ptr<FW::IWriterT<Acts::Surface>>> subWriters;
  std::vector<std::shared_ptr<std::ofstream>>               subStreams;
  // loop and create
  for (auto sdet : subDetectors) {
    // sub detector stream
    auto        sdStream = std::shared_ptr<std::ofstream>(new std::ofstream);
    std::string sdOutputName = sdet + std::string(".obj");
    sdStream->open(sdOutputName);
    // object surface writers
    FWObj::ObjSurfaceWriter::Config sdObjWriterConfig(sdet,
                                                      Acts::Logging::INFO);
    sdObjWriterConfig.filePrefix         = "mtllib materials.mtl";
    sdObjWriterConfig.outputPhiSegemnts  = 72;
    sdObjWriterConfig.outputPrecision    = 6;
    sdObjWriterConfig.outputScalor       = 1.;
    sdObjWriterConfig.outputThickness    = 1.;
    sdObjWriterConfig.outputSensitive    = true;
    sdObjWriterConfig.outputLayerSurface = true;
    sdObjWriterConfig.outputStream       = sdStream;
    auto sdObjWriter
        = std::make_shared<FWObj::ObjSurfaceWriter>(sdObjWriterConfig);
    // call initialize
    sdObjWriter->initialize();
    // push back
    subWriters.push_back(sdObjWriter);
    subStreams.push_back(sdStream);
  }
  // configure the tracking geometry writer
  FWObj::ObjTrackingGeometryWriter::Config tgObjWriterConfig(
      "ObjTrackingGeometryWriter", Acts::Logging::VERBOSE);
  tgObjWriterConfig.surfaceWriters       = subWriters;
  tgObjWriterConfig.filePrefix           = "mtllib materials.mtl";
  tgObjWriterConfig.sensitiveGroupPrefix = "usemtl silicon'\n'";
  tgObjWriterConfig.layerPrefix          = "usemtl support'\n'";
  // the tracking geometry writers
  auto tgObjWriter
      = std::make_shared<FWObj::ObjTrackingGeometryWriter>(tgObjWriterConfig);

  // write the tracking geometry object
  tgObjWriter->write(*(dd4Geometry.get()));

  // --------------------------------------------------------------------------------
  // close the output streams
  for (auto sStreams : subStreams) {
    sStreams->close();
  }
}

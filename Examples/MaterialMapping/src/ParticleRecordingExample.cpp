// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/program_options.hpp>
#include "ACTFW/DD4hepDetector/DD4hepDetectorOptions.hpp"
#include "ACTFW/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/MaterialMapping/OutcomeRecording.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/DD4hepG4/DD4hepToG4Svc.hpp"
#include "ACTFW/Plugins/Root/RootMaterialTrackWriter.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ACTFW/MaterialMapping/OutcomeRecordingOptions.hpp"

#include "ACTFW/Plugins/Geant4/OREventAction.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "ACTFW/Barcode/BarcodeSvc.hpp"

namespace po = boost::program_options;

using PartRec = std::vector<std::vector<FW::Geant4::ParticleRecord>>;
using WriteIt = FW::WriterT<PartRec>;

namespace FW
{
class ParticleRecordWriting : public WriteIt
{
public:

	struct Config
    {
      std::string collection;              ///< particle collection to write
      std::shared_ptr<BarcodeSvc>
             barcodeSvc;          ///< the barcode service to decode (optional)
    };
    
    /// Constructor
    ///
    /// @param cfg Configuration struct
    /// @param level Message level declaration
    ParticleRecordWriting(const Config&        cfg,
                       Acts::Logging::Level level = Acts::Logging::INFO) : WriteIt(cfg.collection, "ParticleRecordWriter", level), m_cfg(cfg) {}

    /// Virtual destructor
    ~ParticleRecordWriting() override {}

    /// End-of-run hook
    ProcessCode
    endRun() final override 
    {
		  return ProcessCode::SUCCESS;
	} // TODO: Maybe write the particles here to file

  protected:
    /// @brief Write method called by the base class
    /// @param [in] context is the algorithm context for event information
    /// @param [in] vertices is the process vertex collection for the
    /// particles to be attached
    ProcessCode
    writeT(const AlgorithmContext&       context,
           const PartRec& vertices) final override
	{
		int eventNr = context.eventNumber;
		
		std::cout << "evtNr: " << eventNr << std::endl;
		for(const auto& evt : vertices)
		{
			std::cout << "New Track: " << std::endl;
			for(const auto& p : evt)
			{
				//~ if(p.volume == "No volume")
					std::cout << "TID: " << p.trackid << ",Parent: " << p.parentid << ", PID: " << p.pdg <<
						", Position: (" << p.position[0] << ", " << p.position[1] << ", " << p.position[2] << 
						"), Momentum: (" << p.momentum[0] << ", " << p.momentum[1] << ", " << p.momentum[2] << "), Volume: " << p.volume << ", Process: " << p.process << std::endl;				
			}
		}
		return ProcessCode::SUCCESS;
	}

  private:
    Config     m_cfg;         ///< The config class
};
}

int
main(int argc, char* argv[])
{
  // Declare the supported program options.
  // Setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addSequencerOptions(desc);
  FW::Options::addOutputOptions(desc);
  FW::Options::addDD4hepOptions(desc);
  FW::Options::addOutcomeRecordingOptions(desc);

  // Parse the options
  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  FW::Sequencer g4sequencer(FW::Options::readSequencerConfig(vm));

  size_t nTracks     = 100;
  int    randomSeed1 = 536235167;
  int    randomSeed2 = 729237523;

  Acts::GeometryContext geoContext;

  // DETECTOR:
  // --------------------------------------------------------------------------------
  // DD4Hep detector definition
  // read the detector config & dd4hep detector
  auto dd4HepDetectorConfig
      = FW::Options::readDD4hepConfig<po::variables_map>(vm);
  auto geometrySvc = std::make_shared<FW::DD4hep::DD4hepGeometryService>(
      dd4HepDetectorConfig);
  std::shared_ptr<const Acts::TrackingGeometry> tGeometry
      = geometrySvc->trackingGeometry(geoContext);

  // DD4Hep to Geant4 conversion
  //
  FW::DD4hepG4::DD4hepToG4Svc::Config dgConfig("DD4hepToG4",
                                               Acts::Logging::INFO);
  dgConfig.dd4hepService = geometrySvc;
  auto dd4hepToG4Svc = std::make_shared<FW::DD4hepG4::DD4hepToG4Svc>(dgConfig);

  // --------------------------------------------------------------------------------
  // Geant4 JOB:
  // --------------------------------------------------------------------------------
  // set up the writer for

  // ---------------------------------------------------------------------------------

  // set up the algorithm writing out the material map
  FW::OutcomeRecording::Config g4rConfig = FW::Options::readOutcomeRecordingConfig<po::variables_map>(vm);
  g4rConfig.geant4Service  = dd4hepToG4Svc;
  g4rConfig.seed1          = randomSeed1;
  g4rConfig.seed2          = randomSeed2;
  g4rConfig.particleCollection = "geant-outcome-tracks";
  
  // create the geant4 algorithm
  auto g4rAlgorithm
      = std::make_shared<FW::OutcomeRecording>(g4rConfig, Acts::Logging::INFO);

  // Output directory
  std::string particleCollection = g4rConfig.particleCollection;

    // Write the propagation steps as ROOT TTree
    FW::ParticleRecordWriting::Config config;
    config.collection  = particleCollection;
    auto writer
        = std::make_shared<FW::ParticleRecordWriting>(
            config);
    g4sequencer.addWriter(writer);

  // Append the algorithm and run
  g4sequencer.addAlgorithm(g4rAlgorithm);
  g4sequencer.run();
}

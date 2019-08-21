// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/program_options.hpp>
#include <memory>

#include "Acts/EventData/TrackParameters.hpp"

#include "ACTFW/EventData/Barcode.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Generators/Pythia8ProcessGenerator.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Options/Pythia8Options.hpp"
#include "ACTFW/Plugins/Csv/CsvParticleWriter.hpp"
#include "ACTFW/Plugins/Root/RootParticleWriter.hpp"
#include "ACTFW/TruthTracking/TrackSelector.hpp"
#include "ACTFW/TruthTracking/TruthVerticesToTracks.hpp"
#include "ACTFW/Utilities/Paths.hpp"

#include "VertexFitAlgorithm.hpp"

using namespace FW;

/// Main vertex fitter example executable
///
/// @param argc The argument count
/// @param argv The argument list
int
main(int argc, char* argv[])
{
  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addPythia8Options(desc);
  Options::addOutputOptions(desc);
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) { return EXIT_FAILURE; }

  auto logLevel = Options::readLogLevel(vm);

  // basic services
  auto rndCfg  = Options::readRandomNumbersConfig(vm);
  auto rnd     = std::make_shared<RandomNumbers>(rndCfg);
  auto barcode = std::make_shared<BarcodeSvc>(BarcodeSvc::Config());

  // Set up event generator producing one single hard collision
  Pythia8Generator::Config hardCfg;
  hardCfg.pdgBeam0  = vm["evg-pdgBeam0"].template as<int>();
  hardCfg.pdgBeam1  = vm["evg-pdgBeam1"].template as<int>();
  hardCfg.cmsEnergy = vm["evg-cmsEnergy"].template as<double>();
  hardCfg.settings  = {vm["evg-hsProcess"].template as<std::string>()};

  auto vtxStdXY = vm["evg-vertex-xy-std"].template as<double>();
  auto vtxStdZ  = vm["evg-vertex-z-std"].template as<double>();

  EventGenerator::Config evgenCfg;
  evgenCfg.generators = {{FixedMultiplicityGenerator{1},
                          GaussianVertexGenerator{vtxStdXY, vtxStdXY, vtxStdZ},
                          Pythia8Generator::makeFunction(hardCfg)}};

  evgenCfg.output        = "generated_particles";
  evgenCfg.randomNumbers = rnd;
  evgenCfg.barcodeSvc    = barcode;

  // Set magnetic field
  Acts::Vector3D bField(0., 0., 2. * Acts::units::_T);

  // Set up TruthVerticesToTracks converter algorithm
  TruthVerticesToTracksAlgorithm::Config trkConvConfig;
  trkConvConfig.randomNumberSvc = rnd;
  trkConvConfig.bField          = bField;
  trkConvConfig.input           = evgenCfg.output;
  trkConvConfig.output          = "all_tracks";

  // Set up track selector
  TrackSelector::Config selectorConfig;
  selectorConfig.input       = trkConvConfig.output;
  selectorConfig.output      = "selected_tracks";
  selectorConfig.absEtaMax   = 2.5;
  selectorConfig.rhoMax      = 4 * Acts::units::_mm;
  selectorConfig.ptMin       = 400. * Acts::units::_MeV;
  selectorConfig.keepNeutral = false;

  // Add the fit algorithm with Billoir fitter
  FWE::VertexFitAlgorithm::Config vertexFitCfg;
  vertexFitCfg.trackCollection = selectorConfig.output;
  vertexFitCfg.bField          = bField;

  Sequencer::Config sequencerCfg = Options::readSequencerConfig(vm);
  Sequencer         sequencer(sequencerCfg);

  sequencer.addReader(std::make_shared<EventGenerator>(evgenCfg, logLevel));

  sequencer.addAlgorithm(std::make_shared<TruthVerticesToTracksAlgorithm>(
      trkConvConfig, logLevel));

  sequencer.addAlgorithm(
      std::make_shared<TrackSelector>(selectorConfig, logLevel));

  sequencer.addAlgorithm(
      std::make_shared<FWE::VertexFitAlgorithm>(vertexFitCfg, logLevel));

  return sequencer.run();
}

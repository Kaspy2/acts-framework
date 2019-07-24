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

#include "ACTFW/Barcode/BarcodeSvc.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Options/Pythia8Options.hpp"
#include "ACTFW/Plugins/Csv/CsvParticleWriter.hpp"
#include "ACTFW/Plugins/Root/RootParticleWriter.hpp"
#include "ACTFW/Random/RandomNumbersOptions.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "ACTFW/Utilities/Paths.hpp"

#include "ACTFW/Generators/ParticleSelector.hpp"
#include "ACTFW/Generators/Pythia8ProcessGenerator.hpp"
#include "ACTFW/TruthTracking/TrackSelector.hpp"
#include "ACTFW/TruthTracking/TruthVerticesToTracks.hpp"
#include "VertexFindingAlgorithm.hpp"

using namespace FW;

/// Main vertex finder example executable
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
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  auto logLevel = Options::readLogLevel(vm);

  // basic services
  RandomNumbersSvc::Config rndCfg  = Options::readRandomNumbersConfig(vm);
  auto                     rnd     = std::make_shared<RandomNumbersSvc>(rndCfg);
  auto                     barcode = std::make_shared<BarcodeSvc>(
      BarcodeSvc::Config(), Acts::getDefaultLogger("BarcodeSvc", logLevel));

  // Set up event generator producing one single hard collision
  EventGenerator::Config evgenCfg = Options::readPythia8Options(vm);
  evgenCfg.output                 = "generated_particles";
  evgenCfg.randomNumbers          = rnd;
  evgenCfg.barcodeSvc             = barcode;

  ParticleSelector::Config ptcSelectorCfg;
  ptcSelectorCfg.input       = evgenCfg.output;
  ptcSelectorCfg.output      = "selected_particles";
  ptcSelectorCfg.absEtaMax   = 2.5;
  ptcSelectorCfg.rhoMax      = 4 * Acts::units::_mm;
  ptcSelectorCfg.ptMin       = 400. * Acts::units::_MeV;
  ptcSelectorCfg.keepNeutral = false;

  // Set magnetic field
  Acts::Vector3D bField(0., 0., 1. * Acts::units::_T);

  // Set up TruthVerticesToTracks converter algorithm
  TruthVerticesToTracksAlgorithm::Config trkConvConfig;
  trkConvConfig.doSmearing      = true;
  trkConvConfig.randomNumberSvc = rnd;
  trkConvConfig.bField          = bField;
  trkConvConfig.input           = ptcSelectorCfg.output;
  trkConvConfig.output          = "all_tracks";

  // Set up track selector
  TrackSelector::Config selectorConfig;
  selectorConfig.input       = trkConvConfig.output;
  selectorConfig.output      = "selectedTracks";
  selectorConfig.absEtaMax   = 2.5;
  selectorConfig.rhoMax      = 4 * Acts::units::_mm;
  selectorConfig.ptMin       = 400. * Acts::units::_MeV;
  selectorConfig.keepNeutral = false;

  // Add the finding algorithm
  FWE::VertexFindingAlgorithm::Config vertexFindingCfg;
  vertexFindingCfg.trackCollection = selectorConfig.output;
  vertexFindingCfg.bField          = bField;

  Sequencer::Config sequencerCfg = Options::readSequencerConfig(vm);
  Sequencer         sequencer(sequencerCfg);

  sequencer.addReader(std::make_shared<EventGenerator>(evgenCfg, logLevel));

  sequencer.addAlgorithm(
      std::make_shared<ParticleSelector>(ptcSelectorCfg, logLevel));

  sequencer.addAlgorithm(std::make_shared<TruthVerticesToTracksAlgorithm>(
      trkConvConfig, logLevel));

  sequencer.addAlgorithm(
      std::make_shared<TrackSelector>(selectorConfig, logLevel));

  sequencer.addAlgorithm(std::make_shared<FWE::VertexFindingAlgorithm>(
      vertexFindingCfg, logLevel));

  return sequencer.run();
}

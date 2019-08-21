// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cstdlib>
#include <memory>

#include <Acts/Utilities/Units.hpp>

#include "ACTFW/EventData/Barcode.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Generators/ParticleSelector.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Options/Pythia8Options.hpp"
#include "ACTFW/Plugins/Csv/CsvParticleWriter.hpp"
#include "ACTFW/Plugins/Root/RootParticleWriter.hpp"
#include "ACTFW/Utilities/Paths.hpp"

using namespace Acts::units;
using namespace FW;

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

  auto              logLevel     = Options::readLogLevel(vm);
  Sequencer::Config sequencerCfg = Options::readSequencerConfig(vm);
  Sequencer         sequencer(sequencerCfg);

  // basic services
  auto rnd
      = std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vm));
  auto barcode = std::make_shared<BarcodeSvc>(BarcodeSvc::Config());

  // event generation w/ process guns
  EventGenerator::Config evgenCfg = Options::readPythia8Options(vm);
  evgenCfg.output                 = "generated_particles";
  evgenCfg.randomNumbers          = rnd;
  evgenCfg.barcodeSvc             = barcode;
  sequencer.addReader(std::make_shared<EventGenerator>(evgenCfg, logLevel));

  // event selection
  ParticleSelector::Config selectorCfg;
  selectorCfg.input  = evgenCfg.output;
  selectorCfg.output = "selected_particles";
  //  selectorCfg.absEtaMax = 2.0;
  //  selectorCfg.ptMin     = 0.5 * _GeV;
  selectorCfg.keepNeutral = false;  // retain only charged particles
  sequencer.addAlgorithm(
      std::make_shared<ParticleSelector>(selectorCfg, logLevel));

  // different output modes
  std::string outputDir = vm["output-dir"].as<std::string>();

  if (vm["output-csv"].as<bool>()) {
    Csv::CsvParticleWriter::Config csvWriterCfg;
    csvWriterCfg.collection     = selectorCfg.output;
    csvWriterCfg.outputDir      = outputDir;
    csvWriterCfg.outputFileName = "particles.csv";
    sequencer.addWriter(
        std::make_shared<Csv::CsvParticleWriter>(csvWriterCfg, logLevel));
  }
  if (vm["output-root"].as<bool>()) {
    Root::RootParticleWriter::Config rootWriterCfg;
    rootWriterCfg.collection = selectorCfg.output;
    rootWriterCfg.filePath   = joinPaths(outputDir, "particles.root");
    rootWriterCfg.barcodeSvc = barcode;
    sequencer.addWriter(
        std::make_shared<Root::RootParticleWriter>(rootWriterCfg, logLevel));
  }

  return sequencer.run();
}

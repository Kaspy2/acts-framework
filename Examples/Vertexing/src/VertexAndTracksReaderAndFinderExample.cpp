// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
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
#include "ACTFW/Plugins/Root/RootVertexAndTracksReader.hpp"
#include "ACTFW/Utilities/Paths.hpp"

#include "VertexFindingAlgorithm.hpp"

using namespace FW;

/// Main executable
///
/// @param argc The argument count
/// @param argv The argument list
int
main(int argc, char* argv[])
{
  using namespace boost::program_options;
  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addOutputOptions(desc);
  desc.add_options()("input",
                     value<std::string>()->default_value(""),
                     "Input root file to read.");
  auto vm = Options::parse(desc, argc, argv);

  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  auto logLevel = Options::readLogLevel(vm);

  // Set file to read data from
  std::string fileString = vm["input"].template as<std::string>();

  if (fileString.empty()) {
    std::cout << "Error: Input file not set." << std::endl;
    return EXIT_FAILURE;
  }

  Root::RootVertexAndTracksReader::Config vtxAndTracksReaderCfg;
  vtxAndTracksReaderCfg.fileList.push_back(fileString);

  // Set magnetic field
  Acts::Vector3D bField(0., 0., 2. * Acts::units::_T);

  // Add the finding algorithm
  FWE::VertexFindingAlgorithm::Config vertexFindingCfg;
  vertexFindingCfg.trackCollection = vtxAndTracksReaderCfg.outputCollection;
  vertexFindingCfg.bField          = bField;

  Sequencer::Config sequencerCfg = Options::readSequencerConfig(vm);
  Sequencer         sequencer(sequencerCfg);

  sequencer.addReader(
      std::make_shared<Root::RootVertexAndTracksReader>(vtxAndTracksReaderCfg));

  sequencer.addAlgorithm(std::make_shared<FWE::VertexFindingAlgorithm>(
      vertexFindingCfg, logLevel));

  return sequencer.run();
}

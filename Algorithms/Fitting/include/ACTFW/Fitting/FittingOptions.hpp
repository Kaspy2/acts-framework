// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include "ACTFW/Utilities/Options.hpp"
#include "FittingAlgorithm.hpp"

namespace po = boost::program_options;

namespace FW {

namespace Options {

  /// @brief read the Fitting options
  ///
  /// Adding Fitting specific options to the Options package
  ///
  /// @tparam aopt_t Type of the options object (API bound to boost)
  ///
  /// @param [in] opt_t The options object where the specific digitization
  /// options are attached to
  template <typename aopt_t>
  void
  addFittingOptions(aopt_t& opt)
  {
    opt.add_options()(
        "fatras-sim-particles",
        po::value<std::string>()->default_value("fatras-particles"),
        "The collection of simulated particles.")(
        "fatras-sim-hits",
        po::value<std::string>()->default_value("fatras-hits"),
        "The collection of simulated hits")(
        "fitted-tracks",
        po::value<std::string>()->default_value("fitted-tracks"),
        "The collection of output tracks");
  }

  /// @brief read the fitter specific options and return a Config file
  ///
  ///@tparam omap_t Type of the options map
  ///@param vm the options map to be read out
  template <typename AMAP, typename kalman_Fitter_t>
  typename FittingAlgorithm<kalman_Fitter_t>::Config
  readFittingConfig(const AMAP& vm, kalman_Fitter_t& fitter)
  {
    // Create a config
    typename FittingAlgorithm<kalman_Fitter_t>::Config fittingConfig(
        std::move(fitter));

    // set the collections
    fittingConfig.simulatedHitCollection
        = vm["fatras-sim-hits"].template as<std::string>();
    fittingConfig.simulatedEventCollection
        = vm["fatras-sim-particles"].template as<std::string>();
    fittingConfig.trackCollection
        = vm["fitted-tracks"].template as<std::string>();

    // and return the config
    return fittingConfig;
  }

}  // namespace Options
}  // namespace FW
